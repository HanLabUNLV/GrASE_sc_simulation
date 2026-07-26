#!/usr/bin/env Rscript
# run_model.R -- thin, model-agnostic driver for the single-cell 3-counts x 3-models
# DTU grid. Reads one count type's `sc` table (gene,event,diff,n,groups; ref = n-diff),
# runs the requested GrASE model + dispersion mode, and writes a results table with
# BH-adjusted p-values and a significance flag. It reuses the package's validated
# scaffold (`run_one_comparison`) and model plugins (`test_model_*`) unchanged; only
# the count builder upstream differs per count type (sj / exon / adjexon).
#
# The three count types all reduce to the same per-event structure, so this one driver
# serves all nine cells of the grid. Dispersion modes:
#   betabinom  : disp EB   (phi_estimate_glmmTMB -> moderate -> test_model_glmmTMB_EB)
#   negbinom   : disp EB   (estimate_theta       -> moderate -> test_model_negbinom)
#   mixedbinom : disp rmats (rMATS logit-normal variance, no estimation) [primary]
#                disp EB    (estimate_logitnormal_var -> moderate -> test_model_mixedbinom)
#
# Significance = BH FDR per contrast (a simple, model-comparable rule; the package's
# nested_BH + independent filtering is bubble-specific and not used here).

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(glmmTMB); library(parallel)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts",  type = "character", help = "sc count table tsv: gene,event,diff,n,groups"),
  make_option("--model",   type = "character", help = "betabinom | negbinom | mixedbinom"),
  make_option("--disp",    type = "character", default = NULL,
              help = "EB | rmats (default: EB, except mixedbinom -> rmats)"),
  make_option("--out",     type = "character", help = "output results tsv"),
  make_option("--pkg_r",   type = "character",
              default = "/mnt/data1/home/mirahan/scGrASE/GrASE/Rpkg/R",
              help = "path to GrASE Rpkg/R (sourced, so the new models are picked up)"),
  make_option("--padj_threshold", type = "double", default = 0.05),
  make_option("--min_effect",     type = "double", default = 0,
              help = "min |effect_size| for the significant flag (0 = off)"),
  make_option("--mc_cores", type = "integer", default = 16L)
)))

stopifnot(!is.null(opt$counts), !is.null(opt$model), !is.null(opt$out))
disp <- opt$disp
if (is.null(disp)) disp <- if (opt$model == "mixedbinom") "rmats" else "EB"

## ---- load the package (source, so the cloned new models are used) ------------
src <- function(f) sys.source(file.path(opt$pkg_r, f), envir = globalenv())
src("exontest_functions.R")      # group_by_event, moderate_phi_log_scale,
                                 # phi_estimate_glmmTMB, test_model_glmmTMB_EB,
                                 # run_one_comparison
src("pvalueAdjustment.R")        # (kept for parity; BH used below)
src("model_negbinom.R")          # estimate_theta, test_model_negbinom
src("model_mixedbinom.R")        # estimate_logitnormal_var, test_model_mixedbinom, rmats var

## ---- read counts -------------------------------------------------------------
sc <- read.table(opt$counts, header = TRUE, sep = "\t",
                 row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
need <- c("gene", "event", "diff", "n", "groups")
if (!all(need %in% names(sc)))
  stop("counts must have columns: ", paste(need, collapse = ", "),
       " (have: ", paste(names(sc), collapse = ", "), ")")
sc <- sc %>%
  filter(!is.na(diff), !is.na(n), n > 0) %>%
  mutate(groups = factor(groups))
# the model plugins read the distinguishing count as dd$y; keep both.
sc$y <- sc$diff
cat(sprintf("[run_model] %s counts: %d rows, %d events, %d groups\n",
            basename(opt$counts), nrow(sc),
            nrow(distinct(sc, gene, event)), nlevels(sc$groups)))

## ---- model dispatch ----------------------------------------------------------
# pick the per-event dispersion estimator (EB only) and the test plugin.
estimator <- switch(opt$model,
  betabinom  = phi_estimate_glmmTMB,
  negbinom   = estimate_theta,
  mixedbinom = estimate_logitnormal_var,
  stop("unknown --model: ", opt$model))
test_fn <- switch(opt$model,
  betabinom  = test_model_glmmTMB_EB,
  negbinom   = test_model_negbinom,
  mixedbinom = test_model_mixedbinom,
  stop("unknown --model: ", opt$model))
model_label <- paste0(opt$model, "_", disp)

# Contrast matrix. betabinom's test requires an explicit L (case-vs-control over the
# `0 + groups` coefficients); negbinom/mixedbinom build their own contrast internally
# when L is NULL (interaction contrast for NB, groups contrast for mixed).
L <- NULL
if (opt$model == "betabinom") {
  levs <- levels(sc$groups)
  if (length(levs) != 2) stop("betabinom driver expects exactly 2 groups")
  L <- matrix(c(-1, 1), nrow = 2,
              dimnames = list(paste0("groups", levs), "group_diff_usage"))
}

## ---- EB moderation (shared moderator) ----------------------------------------
# When disp == EB, estimate the per-event dispersion, moderate on the log scale with
# the SINGLE shared moderator, and join z_mod back onto every row of its event. If the
# counts carry a `comparison` column (adjacent-exon bubbles have two sides per event),
# moderate independently per comparison and key the join on it -- their dispersion
# distributions differ, matching the beta-binomial driver in the package.
has_comp <- "comparison" %in% names(sc)
if (disp == "EB") {
  gd <- group_by_event(sc, "diff", "n")
  cat(sprintf("[run_model] EB: estimating dispersion for %d event-groups...\n", length(gd)))
  phis <- bind_rows(mclapply(gd, function(dd) {
    r <- tryCatch(estimator(dd), error = function(e) NULL)
    if (!is.null(r) && has_comp) r$comparison <- dd$comparison[1]
    r
  }, mc.cores = opt$mc_cores))
  if (is.null(phis) || nrow(phis) == 0) stop("no dispersion estimates succeeded")
  if (has_comp) {
    phi_table <- phis %>% filter(!is.na(comparison)) %>%
      group_by(comparison) %>%
      group_modify(~ moderate_phi_log_scale(.x)) %>% ungroup()
    join_keys <- c("gene", "event", "comparison")
  } else {
    phi_table <- moderate_phi_log_scale(phis)
    join_keys <- c("gene", "event")
  }
  fallback_z <- median(phi_table$z_mod, na.rm = TRUE)
  sc <- sc %>% left_join(phi_table[, c(join_keys, "z_mod")], by = join_keys)
  sc$z_mod[is.na(sc$z_mod)] <- fallback_z
  cat(sprintf("[run_model] EB: %d/%d event-groups moderated (fallback z_mod = %.3f)\n",
              sum(!is.na(phi_table$z_mod)), length(gd), fallback_z))
}

## ---- run the test (reuse the package scaffold) -------------------------------
err_log <- sub("\\.[^.]*$", ".errors.log", opt$out)
extra <- list(sc = sc, test_fn = test_fn, err_log = err_log,
              model_label = model_label, L = L, mc_cores = opt$mc_cores)
if (opt$model == "mixedbinom") extra$disp <- disp     # mixed test takes disp
results <- do.call(run_one_comparison, extra)

if (is.null(results) || nrow(results) == 0) {
  cat("[run_model] WARNING: no results produced\n")
  write.table(data.frame(), file = opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
  quit(save = "no", status = 0)
}

## ---- FDR + significance flag -------------------------------------------------
results <- results %>%
  group_by(contrast) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  ungroup()
results$significant <- results$padj < opt$padj_threshold &
  (opt$min_effect == 0 | abs(results$effect_size) >= opt$min_effect)

write.table(results, file = opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[run_model] wrote %d results (%d significant) -> %s\n",
            nrow(results), sum(results$significant, na.rm = TRUE), opt$out))
