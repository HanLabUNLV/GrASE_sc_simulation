# NB normalization test: effect of scRNA-seq size-factor normalization

Question: does adding scRNA-seq library-size normalization (a size-factor offset) to the
NB DTU model change its performance?

Normalization = library-size size factors (sf_cell = total exonic counts / mean), the
practical scRNA-seq standard for NB offsets (glmGamPoi / Seurat default; NOT DESeq2
median-of-ratios, which breaks on scRNA-seq zeros). Computed from the per-cell DEXSeq
total counts (range 3.86M-9.12M, sf 0.644-1.522). Exon arm, ALL sim_type, restricted,
padj 0.05.

| norm mode    | formula                              | sens  | spec  | prec  | TP   | FP   | fit speed          |
|--------------|--------------------------------------|-------|-------|-------|------|------|--------------------|
| none         | count ~ 0+groups:part + (1|id)       | 0.347 | 0.990 | 0.498 | 1653 | 1666 | slow (integration) |
| offset_id    | ... + (1|id) + offset(log sf)        | 0.347 | 0.990 | 0.499 | 1653 | 1657 | slow               |
| offset_only  | ... + offset(log sf)   (no RE)       | 0.136 | 0.997 | 0.549 | 647  | 532  | fast (no integration) |

## Findings

1. Size-factor normalization ON TOP of the (1|id) cell effect has ZERO effect
   (0.347 -> 0.347, TP 1653 = 1653). The random effect already absorbs per-cell scale, so
   an explicit library-size offset is fully redundant. Makes sense: for the DTU
   groups:part INTERACTION, cell-level scaling cancels.
2. Replacing (1|id) with the offset (offset_only) HALVES sensitivity (0.347 -> 0.136),
   with a small specificity/precision gain. The random effect does essential work beyond
   normalization -- it captures cell-level overdispersion and the diff/ref pairing the DTU
   test relies on; a fixed offset supplies scale but not that.
3. offset_only is the no-random-effect / no-Laplace variant (~5 min estimation vs ~25 min)
   -- confirming the integration can be removed, but at ~60% of sensitivity. Not worth it.

## Conclusion

For scRNA-seq NB DTU, explicit size-factor normalization adds nothing given the cell
random effect, and the random effect cannot be cheaply swapped for a size-factor offset
without losing half the sensitivity. Keep the (1|id) default; no normalization needed.
(A rigorous follow-up would swap library-size for scran pooled size factors, but since the
offset is redundant on top of (1|id), the choice of size factor is unlikely to matter.)

How it was run: benchmark/stage7f_nb_norm.sh; run_model.R --norm {offset_id,offset_only}
--sizefactors work/grid/counts/exon_sizefactors.txt; model_negbinom.R .nb_formula.
