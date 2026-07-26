#!/usr/bin/env Rscript
# build_sj_counts.R -- SJ (splice-junction) count builder for the grid's rMATS-proxy
# arm. Turns rMATS-turbo's per-cell inclusion/skip junction counts (*.MATS.JCEC.txt)
# into the shared sc table (gene, event, diff, n, groups) that run_model.R consumes,
# so the SJ counts are tested by the SAME glmmTMB models as the exon and adjacent-exon
# counts. This is count-building only (adapted from jaquino's rmats_cnts.R); the model
# and the FDR live in run_model.R.
#
# Per event, rMATS stores one comma-delimited value per cell:
#   IJC_SAMPLE_1 / SJC_SAMPLE_1  inclusion / skip counts, group1 (control)
#   IJC_SAMPLE_2 / SJC_SAMPLE_2  inclusion / skip counts, group2 (case)
# We emit, per (event x cell): diff = inclusion, n = inclusion + skip.
#
# event id = "<event_type>_<ID>" and gene = GeneID, so the SJ eval can split the id
# back and join to the junction GT (sim_junction_gt.txt) by (event_type, ID) -- the
# same key the rMATS-native tool eval uses.

suppressPackageStartupMessages({ library(optparse); library(dplyr) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--rmats_dir", type = "character",
              help = "rMATS post dir with *.MATS.JCEC.txt (rmats_post_group1_group2)"),
  make_option("--out", type = "character", help = "output sc counts tsv"),
  make_option("--types", type = "character", default = "SE,A3SS,A5SS,RI",
              help = "event types to include [default: %default] (MXE excluded: dual inclusion form)"),
  make_option("--min_n", type = "integer", default = 1L,
              help = "drop (event x cell) rows with total n below this [default: 1]")
)))
stopifnot(!is.null(opt$rmats_dir), !is.null(opt$out))
types <- strsplit(opt$types, ",")[[1]]

# split one event's comma-delimited per-cell count string into an integer vector.
csv_int <- function(x) as.integer(strsplit(as.character(x), ",", fixed = TRUE)[[1]])

# build the long per-cell rows for one group of one event type.
long_group <- function(df, ijc_col, sjc_col, grp_label) {
  ijc <- strsplit(as.character(df[[ijc_col]]), ",", fixed = TRUE)
  sjc <- strsplit(as.character(df[[sjc_col]]), ",", fixed = TRUE)
  ncells <- lengths(ijc)
  inc <- as.integer(unlist(ijc))
  skp <- as.integer(unlist(sjc))
  data.frame(
    gene   = rep(df$gene, ncells),
    event  = rep(df$event, ncells),
    diff   = inc,
    n      = inc + skp,
    groups = grp_label,
    stringsAsFactors = FALSE)
}

all_rows <- list()
for (etype in types) {
  f <- file.path(opt$rmats_dir, paste0(etype, ".MATS.JCEC.txt"))
  if (!file.exists(f)) { cat(sprintf("[build_sj] skip %s (no file)\n", etype)); next }
  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  names(df) <- make.unique(names(df))
  df$GeneID <- gsub('"', '', df$GeneID)
  df$gene  <- df$GeneID
  df$event <- paste0(etype, "_", df$ID)
  g1 <- long_group(df, "IJC_SAMPLE_1", "SJC_SAMPLE_1", "group1")   # control
  g2 <- long_group(df, "IJC_SAMPLE_2", "SJC_SAMPLE_2", "group2")   # case
  all_rows[[etype]] <- bind_rows(g1, g2)
  cat(sprintf("[build_sj] %s: %d events -> %d (event x cell) rows\n",
              etype, nrow(df), nrow(all_rows[[etype]])))
}
sc <- bind_rows(all_rows)
sc <- sc[!is.na(sc$n) & sc$n >= opt$min_n & !is.na(sc$diff) & sc$diff >= 0, ]
# groups control-first so downstream contrasts read case - control.
sc$groups <- factor(sc$groups, levels = c("group1", "group2"))

write.table(sc, file = opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[build_sj] wrote %d rows, %d events -> %s\n",
            nrow(sc), nrow(distinct(sc, gene, event)), opt$out))
