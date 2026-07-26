#!/usr/bin/env Rscript
# NB design audit: jaquino's exact adjexonnegbinomLRT (fixed id + LRT + free theta) vs
# mine (random (1|id) + Wald + free theta) on the SAME adjacent-exon events. Isolate the
# design effect on detection, using betabinom-significant events as a high-confidence
# reference. Free theta on both to remove the EB confound; compares id+test design only.
suppressPackageStartupMessages({ library(dplyr); library(glmmTMB) })
W <- "/mnt/data1/home/mirahan/scGrASE"
src <- file.path(W, "GrASE/Rpkg/R")
sys.source(file.path(src, "exontest_functions.R"), envir = globalenv())

# internal bipartition counts -> d1r comparison (this=diff1, others=ref)
raw <- read.table(file.path(W, "work/grase/bipartition.internal.counts/bipartition.internal.exoncnt.combined.txt"),
                  header = TRUE, row.names = NULL, check.names = FALSE)
sc <- raw %>% filter(!is.na(diff1) & !is.na(ref)) %>%
  transmute(gene, event, groups = factor(groups), this = diff1, others = ref)
gd <- sc %>% group_by(gene, event) %>% group_split()

# betabinom-significant events (high-confidence reference)
bb <- read.table(file.path(W, "work/grase/bipartition.test/test_bipartition.internal_betabinom_EBapprox.annotated.txt"),
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)
bb_sig <- unique(paste(bb$gene, bb$event)[bb$significant %in% c(TRUE, "TRUE")])

# jaquino's exact model (fixed id + LRT, free theta)
jaq_nb <- function(d) {
  d <- d[(d$this + d$others) > 0, ]
  if (nrow(d) < 4 || length(unique(d$groups)) < 2) return(NA_real_)
  n <- nrow(d)
  dat <- data.frame(y = c(d$this, d$others),
                    id = factor(rep(seq_len(n), 2)),
                    thot = factor(rep(c("this","others"), each = n), levels = c("others","this")),
                    groups = factor(rep(as.character(d$groups), 2)))
  ff <- tryCatch(glmmTMB(y ~ id + thot + groups:thot, data = dat, family = nbinom2, se = FALSE), error = function(e) NULL)
  fr <- tryCatch(glmmTMB(y ~ id + thot,               data = dat, family = nbinom2, se = FALSE), error = function(e) NULL)
  if (is.null(ff) || is.null(fr) || is.na(logLik(ff)) || is.na(logLik(fr))) return(NA_real_)
  a <- tryCatch(anova(ff, fr, test = "LRT"), error = function(e) NULL)
  if (is.null(a)) return(NA_real_)
  a$`Pr(>Chisq)`[2]
}

# my model (random (1|id) + Wald, free theta) -- same structure as test_model_negbinom
# but theta free (no EB) to isolate the id/test design
mine_nb <- function(d) {
  d <- d[(d$this + d$others) > 0, ]
  if (nrow(d) < 4 || length(unique(d$groups)) < 2) return(NA_real_)
  n <- nrow(d)
  ld <- data.frame(count = c(d$this, d$others),
                   part = factor(rep(c("diff","ref"), each = n), levels = c("ref","diff")),
                   groups = factor(rep(as.character(d$groups), 2)),
                   id = factor(rep(seq_len(n), 2)))
  ld$groups <- droplevels(ld$groups)
  m <- tryCatch(glmmTMB(count ~ 0 + groups:part + (1|id), data = ld, family = nbinom2, se = TRUE), error = function(e) NULL)
  if (is.null(m) || is.na(logLik(m))) return(NA_real_)
  b <- glmmTMB::fixef(m)$cond; V <- vcov(m)$cond; cn <- names(b)
  dd <- grep(":partdiff$", cn, value = TRUE); rr <- grep(":partref$", cn, value = TRUE)
  if (length(dd) != 2 || length(rr) != 2) return(NA_real_)
  L <- setNames(numeric(length(cn)), cn)
  L[dd[2]] <- 1; L[rr[2]] <- -1; L[dd[1]] <- -1; L[rr[1]] <- 1
  est <- sum(L * b); se <- sqrt(as.numeric(t(L) %*% V %*% L))
  if (!is.finite(se) || se <= 0) return(NA_real_)
  pchisq((est/se)^2, df = 1, lower.tail = FALSE)
}

set.seed <- NULL
# sample: all betabinom-sig events present + a random 200 others (bounded run time)
keys <- sapply(gd, function(d) paste(d$gene[1], d$event[1]))
is_sig <- keys %in% bb_sig
idx <- c(which(is_sig), head(which(!is_sig), 200))
idx <- idx[seq_len(min(length(idx), 260))]
cat("auditing", length(idx), "events (", sum(keys[idx] %in% bb_sig), "betabinom-sig )\n")

res <- lapply(idx, function(i) {
  d <- gd[[i]]
  data.frame(key = keys[i], bb_sig = keys[i] %in% bb_sig,
             p_jaq = jaq_nb(d), p_mine = mine_nb(d), stringsAsFactors = FALSE)
})
r <- bind_rows(res)
r <- r[!is.na(r$p_jaq) & !is.na(r$p_mine), ]
cat("\nfit on", nrow(r), "events\n")
cat(sprintf("raw p<0.05:  jaquino(fixed id+LRT)=%d   mine(random id+Wald)=%d\n",
            sum(r$p_jaq < 0.05), sum(r$p_mine < 0.05)))
cat("\namong betabinom-significant events (high-confidence true):\n")
s <- r[r$bb_sig, ]
cat(sprintf("  n=%d;  jaquino detects %d (%.0f%%);  mine detects %d (%.0f%%)\n",
            nrow(s), sum(s$p_jaq < 0.05), 100*mean(s$p_jaq < 0.05),
            sum(s$p_mine < 0.05), 100*mean(s$p_mine < 0.05)))
cat("\nmedian p on betabinom-sig events:  jaquino=%.2g  mine=%.2g\n")
cat(sprintf("  jaquino median p=%.3g  mine median p=%.3g\n", median(s$p_jaq), median(s$p_mine)))
