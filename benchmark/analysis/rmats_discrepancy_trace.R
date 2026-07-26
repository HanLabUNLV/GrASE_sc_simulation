#!/usr/bin/env Rscript
# Trace one junction where the proxy calls DTU but native rMATS does not (GT-positive),
# dumping both models' inputs and outputs so we can see where the test diverges.
suppressPackageStartupMessages({ library(dplyr) })
W <- "/mnt/data1/home/mirahan/scGrASE"

proxy <- read.table(file.path(W,"work/grid/results/sj__mixedbinom_rmats.txt"),
                    header=TRUE, sep="\t", stringsAsFactors=FALSE)
proxy$event_type <- sub("_.*$","",proxy$event); proxy$ID <- sub("^[^_]*_","",proxy$event)
gt <- read.table(file.path(W,"work/truth/sim_junction_gt_matched.txt"),
                 header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
gt$ID <- as.character(gt$ID)
se <- read.table(file.path(W,"work/rMATS/rmats_post_group1_group2/SE.MATS.JCEC.txt"),
                 header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
se$GeneID <- gsub('"','',se$GeneID); se$ID <- as.character(se$ID)

# candidates: SE, gt_positive, proxy padj<0.01, native FDR>0.5
p <- proxy %>% filter(event_type=="SE") %>% mutate(ID=as.character(ID))
cand <- p %>% inner_join(gt %>% filter(event_type=="SE", gt_positive) %>% select(ID,sim_type,dPSI), by="ID") %>%
              inner_join(se %>% select(ID,FDR,PValue,IncLevel1,IncLevel2,IncLevelDifference,
                                       IJC_SAMPLE_1,SJC_SAMPLE_1,IJC_SAMPLE_2,SJC_SAMPLE_2), by="ID") %>%
              filter(padj < 0.01, FDR > 0.5) %>% arrange(padj)
cat("discrepancy candidates (proxy sig, rMATS not, GT+):", nrow(cand), "\n\n")
if (nrow(cand)==0) quit(save="no")

e <- cand[1,]
csv <- function(x) as.integer(strsplit(as.character(x),",",fixed=TRUE)[[1]])
ijc1<-csv(e$IJC_SAMPLE_1); sjc1<-csv(e$SJC_SAMPLE_1); ijc2<-csv(e$IJC_SAMPLE_2); sjc2<-csv(e$SJC_SAMPLE_2)
psi1 <- ijc1/(ijc1+sjc1); psi2 <- ijc2/(ijc2+sjc2)

cat("=== event", e$event, "gene", proxy$gene[match(e$event,proxy$event)], "===\n")
cat(sprintf("GT: gt_positive=TRUE sim_type=%s dPSI(gt)=%.3f\n", e$sim_type, e$dPSI))
cat(sprintf("\n-- counts (per cell) --\n g1 inclusion mean=%.1f skip mean=%.1f  n cells=%d\n",
            mean(ijc1), mean(sjc1), length(ijc1)))
cat(sprintf(" g2 inclusion mean=%.1f skip mean=%.1f  n cells=%d\n", mean(ijc2), mean(sjc2), length(ijc2)))
cat(sprintf(" mean PSI: g1=%.3f g2=%.3f  observed dPSI=%.3f\n",
            mean(psi1,na.rm=TRUE), mean(psi2,na.rm=TRUE), mean(psi2,na.rm=TRUE)-mean(psi1,na.rm=TRUE)))

cat("\n-- NATIVE rMATS --\n")
cat(sprintf(" IncLevel1=%.3f IncLevel2=%.3f IncLevelDifference=%.3f\n",
            mean(csv(gsub('NA','',e$IncLevel1))[!is.na(csv(gsub('NA','',e$IncLevel1)))]),
            mean(csv(gsub('NA','',e$IncLevel2))[!is.na(csv(gsub('NA','',e$IncLevel2)))]), e$IncLevelDifference))
cat(sprintf(" PValue=%.3g  FDR=%.3g  -> NOT significant\n", e$PValue, e$FDR))

cat("\n-- PROXY (mixedbinom_rmats) --\n")
cat(sprintf(" p.value=%.3g padj=%.3g LRT(z^2)=%.2f phi(var used)=%.4f effect(logit)=%.3f -> SIGNIFICANT\n",
            e$p.value, e$padj, e$LRT, e$phi, e$effect_size))
# recompute the proxy's rMATS variance for this event
psi_all <- c(psi1,psi2); g <- c(rep("g1",length(psi1)),rep("g2",length(psi2)))
psi_all <- psi_all[is.finite(psi_all)]; g <- g[is.finite(c(psi1,psi2))]
ss<-0; df<-0; for(lv in unique(g)){pp<-psi_all[g==lv]; if(length(pp)>=2){ss<-ss+sum((pp-mean(pp))^2); df<-df+length(pp)-1}}
var_emp <- ss/df; var_reg <- max(0.01, 10*var_emp)
cat(sprintf(" empirical within-group var(psi)=%.4f  -> rMATS reg var=max(0.01,10x)=%.4f\n", var_emp, var_reg))
cat(sprintf(" (proxy phi reported=%.4f)\n", e$phi))
