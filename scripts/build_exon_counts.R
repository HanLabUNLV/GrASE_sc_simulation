#!/usr/bin/env Rscript
# build_exon_counts.R -- exon-bin count builder for the grid's DEXSeq-proxy arm.
# Turns DEXSeq per-cell HTSeq bin counts into the shared sc table (gene, event, diff,
# n, groups) so the exon-bin counts run through the SAME glmmTMB models as the SJ and
# adjacent-exon counts. This IS the scalable DEXSeq arm: it builds the DTU counts
# directly (bin vs rest-of-gene) and hands them to run_model.R, avoiding classic
# testForDEU, which does not scale to 160 samples. Count-building only (the DTU
# structure of jaquino's dexseq_sim.R: y = bin, n = y + other-bins-in-gene = gene total).
#
# Input: DEXSeq count files, one per cell, each "GENE:bin<TAB>count" (dexseq_count.py).
#   group1 dir = control cells, group2 dir = case cells.
# Per (bin x cell): diff = bin count, n = the cell's total over the bin's gene.
# event = the full "GENE:bin" feature id, gene = the GENE part.
#
# Genes never simulated (zero across all cells) and all-zero bins are dropped: they are
# structurally untestable and would only bloat the table; the eval handles the universe.

suppressPackageStartupMessages({ library(optparse); library(dplyr) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--group1_dir", type = "character", help = "control count files dir"),
  make_option("--group2_dir", type = "character", help = "case count files dir"),
  make_option("--out", type = "character", help = "output sc counts tsv"),
  make_option("--min_gene_total", type = "integer", default = 1L,
              help = "drop (bin x cell) rows whose gene total n is below this [default: 1]")
)))
stopifnot(!is.null(opt$group1_dir), !is.null(opt$group2_dir), !is.null(opt$out))

f1 <- sort(list.files(opt$group1_dir, pattern = "\\.txt$", full.names = TRUE))
f2 <- sort(list.files(opt$group2_dir, pattern = "\\.txt$", full.names = TRUE))
files  <- c(f1, f2)
groups <- c(rep("group1", length(f1)), rep("group2", length(f2)))   # control, case
cat(sprintf("[build_exon] %d control + %d case count files\n", length(f1), length(f2)))

# read one count file -> named integer vector (feature -> count), dropping the
# HTSeq summary rows (_ambiguous / _empty / _lowaqual / _notaligned ...).
read_counts <- function(path) {
  d <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                  col.names = c("feature", "count"), quote = "")
  d <- d[!grepl("^_", d$feature), ]
  setNames(as.integer(d$count), d$feature)
}

# master feature list from the first file; align every file to it by name so a
# differing row order never mis-pairs counts.
master <- names(read_counts(files[1]))
M <- matrix(0L, nrow = length(master), ncol = length(files),
            dimnames = list(master, NULL))
for (j in seq_along(files)) {
  v <- read_counts(files[j])
  M[, j] <- v[master]                       # NA if a feature is absent -> set to 0
  if (j %% 40 == 0) cat(sprintf("[build_exon] read %d/%d files\n", j, length(files)))
}
M[is.na(M)] <- 0L

gene <- sub(":.*$", "", master)
# per-cell gene totals, then map each bin to its gene's total (n = gene total).
gene_tot <- rowsum(M, group = gene)          # (ngenes x ncells)
Gn <- gene_tot[gene, , drop = FALSE]         # (nbins x ncells), n for each bin

# drop structurally untestable bins: all-zero bin, or gene never expressed.
keep <- rowSums(M) > 0 & rowSums(Gn) > 0
cat(sprintf("[build_exon] %d bins total, %d kept (nonzero bin & expressed gene)\n",
            length(master), sum(keep)))
M <- M[keep, , drop = FALSE]; Gn <- Gn[keep, , drop = FALSE]
feat <- master[keep]; gene <- gene[keep]

ncells <- length(files)
sc <- data.frame(
  gene   = rep(gene, ncells),
  event  = rep(feat, ncells),
  diff   = as.integer(M),
  n      = as.integer(Gn),
  groups = rep(groups, each = nrow(M)),
  stringsAsFactors = FALSE)
sc <- sc[!is.na(sc$n) & sc$n >= opt$min_gene_total & !is.na(sc$diff) & sc$diff >= 0, ]
sc$groups <- factor(sc$groups, levels = c("group1", "group2"))   # control first

write.table(sc, file = opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[build_exon] wrote %d rows, %d bins (events) -> %s\n",
            nrow(sc), nrow(distinct(sc, gene, event)), opt$out))
