#!/usr/bin/env Rscript
# build_adjexon_counts.R -- adjacent-exon count builder for the grid's GrASE arm.
# Turns GrASE's bubble bipartition counts (exoncnt.combined.txt: diff1, diff2, ref per
# cell) into the shared sc table so the adjacent-exon counts run through the SAME
# glmmTMB models as the SJ and exon counts. This replicates exontest.R's
# `make_comparison`: each bubble contributes two comparisons -- diff1_vs_ref and
# diff2_vs_ref -- stacked with a `comparison` column, which group_by_event/run_model.R
# key on so the two sides stay separate. diff = the distinguishing path count, n = diff
# + ref. Multiple bipartition files (internal, TSSTTS) may be combined.

suppressPackageStartupMessages({ library(optparse); library(dplyr); library(data.table) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts", type = "character",
              help = "comma-separated bipartition exoncnt.combined.txt files (internal[,TSSTTS])"),
  make_option("--out", type = "character", help = "output sc counts tsv"),
  make_option("--min_n", type = "integer", default = 1L,
              help = "drop (event x cell) rows with total n below this [default: 1]")
)))
stopifnot(!is.null(opt$counts), !is.null(opt$out))
files <- strsplit(opt$counts, ",")[[1]]

make_comparison <- function(sc, diff_col, comp_name) {
  sc %>%
    filter(!is.na(.data[[diff_col]]) & !is.na(ref)) %>%
    transmute(gene, event, groups,
              diff       = .data[[diff_col]],
              n          = .data[[diff_col]] + ref,
              comparison = comp_name)
}

all_rows <- list()
for (f in files) {
  if (!file.exists(f)) { cat(sprintf("[build_adjexon] skip (missing): %s\n", f)); next }
  # fread: the bipartition tables can be tens of millions of rows (TSSTTS). The files
  # have an unnamed leading row-index column (13 fields vs 12 header names), so fread
  # auto-adds a first name; read all, then select the columns we need by name.
  raw <- as.data.frame(data.table::fread(f, sep = "\t"))   # auto-header handles the shift
  need <- c("gene", "event", "groups", "ref", "diff1", "diff2")
  if (!all(need %in% names(raw)))
    stop("bipartition file missing columns (", basename(f), "): need ",
         paste(need, collapse = ", "))
  d1r <- make_comparison(raw, "diff1", "diff1_vs_ref")
  d2r <- make_comparison(raw, "diff2", "diff2_vs_ref")
  all_rows[[f]] <- bind_rows(d1r, d2r)
  cat(sprintf("[build_adjexon] %s: %d rows (%d events)\n",
              basename(f), nrow(all_rows[[f]]), nrow(distinct(all_rows[[f]], gene, event))))
}
sc <- bind_rows(all_rows)
sc <- sc[!is.na(sc$n) & sc$n >= opt$min_n & !is.na(sc$diff) & sc$diff >= 0, ]
sc$groups <- factor(sc$groups)

write.table(sc, file = opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[build_adjexon] wrote %d rows, %d events (2 comparisons/bubble) -> %s\n",
            nrow(sc), nrow(distinct(sc, gene, event)), opt$out))
