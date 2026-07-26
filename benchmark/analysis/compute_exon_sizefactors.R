g1 <- sort(list.files("work/DEXSeq/count_files/group1", pattern="\\.txt$", full.names=TRUE))
g2 <- sort(list.files("work/DEXSeq/count_files/group2", pattern="\\.txt$", full.names=TRUE))
files <- c(g1, g2)
totals <- sapply(files, function(f){
  d <- read.table(f, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  sum(as.numeric(d$V2[!grepl("^_", d$V1)]))   # total exonic counts, drop HTSeq summary rows
})
sf <- totals / mean(totals)                    # library-size size factors, mean 1
out <- data.frame(cell = seq_along(sf), library_size = as.integer(totals), sf = round(sf,5))
write.table(out, "work/grid/counts/exon_sizefactors.txt", sep="\t", quote=FALSE, row.names=FALSE)
cat("cells:", length(sf), " lib size range:", min(totals), "-", max(totals),
    " sf range:", round(min(sf),3), "-", round(max(sf),3), "\n")
