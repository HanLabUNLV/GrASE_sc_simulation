#!/usr/bin/env Rscript
# evaluate_grid_sj.R -- evaluate the grid's SJ (splice-junction) arm against the
# junction ground truth. The SJ counts run through the three glmmTMB models are the
# design-intended rMATS proxy; this scores each model's calls vs sim_junction_gt.txt
# and reports sens/spec/precision/recall so they can be cross-checked against jaquino's
# diff_params_4 SJ rows.
#
# Grid SJ results carry event = "<event_type>_<ID>" (from build_sj_counts.R); we split
# that back and join to the junction GT by (event_type, ID) -- the same key the
# rMATS-native tool eval uses. An event is CALLED if its (min over rows) padj < thr.
#
# Universes:
#   full       : every junction-GT event; results absent from the grid are forced
#                non-significant (GT+ -> FN, GT- -> TN), so nothing is silently dropped.
#   restricted : only events actually tested by the grid (present in the results).
#
# Usage: evaluate_grid_sj.R --results f1[,f2,...] --junction_gt sim_junction_gt.txt \
#          --out summary.txt [--padj 0.05] [--dtu_only]

suppressPackageStartupMessages({ library(optparse); library(dplyr); library(tidyr) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--results", type = "character",
              help = "comma-separated grid SJ result files (sj__<model>_<disp>.txt)"),
  make_option("--junction_gt", type = "character", help = "sim_junction_gt.txt"),
  make_option("--out", type = "character", help = "output summary tsv"),
  make_option("--padj", type = "double", default = 0.05),
  make_option("--dtu_only", action = "store_true", default = FALSE,
              help = "restrict FP accounting to DTU-simtype genes (like jaquino's DTU rows)")
)))
stopifnot(!is.null(opt$results), !is.null(opt$junction_gt), !is.null(opt$out))

jgt <- read.table(opt$junction_gt, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "")
jgt$ID         <- as.character(jgt$ID)
jgt$event_type <- as.character(jgt$event_type)
jgt$gt_positive <- as.logical(jgt$gt_positive)
if (opt$dtu_only && "sim_type" %in% names(jgt))
  jgt <- jgt[jgt$sim_type == "DTU", ]
cat(sprintf("[eval_sj] junction GT: %d events, %d gt_positive%s\n",
            nrow(jgt), sum(jgt$gt_positive, na.rm = TRUE),
            if (opt$dtu_only) " (DTU genes only)" else ""))

score_one <- function(res_file) {
  model <- sub("^sj__", "", sub("\\.txt$", "", basename(res_file)))
  r <- read.table(res_file, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "")
  if (nrow(r) == 0) return(NULL)
  # min padj per event (an event is called if any contrast/comparison row is sig).
  r <- r %>%
    separate(event, into = c("event_type", "ID"), sep = "_", extra = "merge") %>%
    group_by(event_type, ID) %>%
    summarise(padj = min(padj, na.rm = TRUE), .groups = "drop")
  r$ID <- as.character(r$ID)

  m <- merge(jgt[, c("event_type", "ID", "gt_positive")], r,
             by = c("event_type", "ID"), all.x = TRUE)
  m$tested <- !is.na(m$padj)
  m$padj[is.na(m$padj)] <- 1               # untested -> forced non-significant
  m$called <- m$padj < opt$padj

  tab <- function(sel) {
    s <- m[sel, ]
    TP <- sum(s$called & s$gt_positive)
    FP <- sum(s$called & !s$gt_positive)
    FN <- sum(!s$called & s$gt_positive)
    TN <- sum(!s$called & !s$gt_positive)
    data.frame(
      TP = TP, FP = FP, FN = FN, TN = TN, total = nrow(s),
      sensitivity = TP / (TP + FN),
      specificity = TN / (TN + FP),
      precision   = if (TP + FP > 0) TP / (TP + FP) else NA_real_,
      recall      = TP / (TP + FN))
  }
  bind_rows(
    cbind(model = model, universe = "full",       tab(rep(TRUE, nrow(m)))),
    cbind(model = model, universe = "restricted", tab(m$tested)))
}

res_files <- strsplit(opt$results, ",")[[1]]
summ <- bind_rows(lapply(res_files, score_one))
summ <- summ[order(summ$universe, summ$model), ]
write.table(summ, file = opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("[eval_sj] wrote", opt$out, "\n\n")
print(summ, row.names = FALSE, digits = 3)
