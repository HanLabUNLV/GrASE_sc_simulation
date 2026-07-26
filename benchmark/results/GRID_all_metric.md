# scGrASE grid -- ALL metric (nulls + DGE charged as FP), confirmed pipeline (10M/80, 50bp)

Every cell holds the confirmed evaluation constant and varies only the model:
- adjacent-exon: exontest.R (new --model branches) -> evaluate_bipartition_test.R vs
  sim_exon_info exon-part GT.
- exon: run_model.R (flat bins, confirmed modules + BH) -> confirmed dexseq_file shape ->
  evaluate_tools_dexseq_vs_gt.R DEXSeq arm vs sim_exon_info exon-part GT.
Validation: adjacent-exon betabinom reproduces the confirmed 0.785/0.216 (DTU) exactly.

sim_type = ALL charges false positives from null (Background) and DGE genes, so precision
is the honest whole-experiment precision. Restricted universe, padj 0.05.

## ALL metric, restricted, @0.05

| count arm     | model             | precision | recall | TP   | FP    | FN   |
|---------------|-------------------|-----------|--------|------|-------|------|
| adjacent-exon | betabinom_EB      | 0.754     | 0.216  | 150  | 49    | 544  |
| adjacent-exon | negbinom_EB       | 0.689     | 0.193  | 133  | 60    | 558  |
| adjacent-exon | mixedbinom_EB     | 0.759     | 0.205  | 142  | 45    | 550  |
| adjacent-exon | mixedbinom_rmats  | 0.232     | 0.444  | 308  | 1022  | 386  |
| exon          | betabinom_EB      | 0.574     | 0.389  | 1979 | 1471  | 3107 |
| exon          | negbinom_EB       | 0.498     | 0.347  | 1653 | 1666  | 3115 |
| exon          | mixedbinom_EB     | 0.560     | 0.357  | 1702 | 1340  | 3066 |
| exon          | mixedbinom_rmats  | 0.132     | 0.781  | 3974 | 26175 | 1112 |

## Sensitivity/specificity vs jaquino thesis (10M/80, @0.05)

sens = recall; spec = TN/(TN+FP). Restricted universe, ALL sim_type.

| arm      | model            | my sens | thesis sens | my spec | thesis spec |
|----------|------------------|---------|-------------|---------|-------------|
| adj-exon | betabinom        | 0.216   | 0.270       | 0.997   | 0.914       |
| adj-exon | mixedbinom_EB    | 0.205   | 0.245       | 0.997   | 0.942       |
| adj-exon | negbinom         | 0.192   | 0.346       | 0.997   | 0.834       |
| adj-exon | mixedbinom_rmats | 0.444   | (new)       | 0.942   | (new)       |
| exon     | betabinom        | 0.415   | 0.625 (DEXSeq all-bins) | 0.991 | 0.687 |
| exon     | negbinom         | 0.347   | -           | 0.990   | -           |
| exon     | mixedbinom_EB    | 0.357   | -           | 0.992   | -           |
| exon     | mixedbinom_rmats | 0.833   | (new)       | 0.838   | (new)       |

Comparable in regime; sensitivity same order, generally a bit lower (betabinom/mixedbinom
map closely). Two caveats: (1) NEGBINOM diverges most (0.192 vs 0.346 adj) -- likely the
random (1|id) vs jaquino's fixed-id design change (audit item). (2) Specificity is NOT
apples-to-apples: my restricted-universe ALL has a large null TN pool so spec sits ~0.99
vs the thesis 0.83-0.94; a different universe/denominator, not a real model difference.
mixedbinom_rmats is a new rMATS-faithful variant not in the thesis.

## Reading it

- The three EB models (beta / nb / mixed) are the high-precision / modest-recall regime;
  mixedbinom_rmats trades precision for recall and, under ALL (nulls charged), its
  precision collapses (adjacent-exon 0.23, exon 0.13 with 26k FP). The rMATS-variance
  floor (0.01) is what drives that.
- Exon (DEXSeq-bin) counts are more sensitive than adjacent-exon (betabinom recall
  0.39 vs 0.22), at lower precision -- consistent with jaquino's thesis (DEXSeq exon bins
  more sensitive than adjacent-exon).
- Native rMATS is not a valid comparator here (degenerate: 99% PValue=1 from its
  zero-count filter on sparse single-cell data -- see RMATS_proxy_vs_native.md).

## How it was run
benchmark/stage7b/7c (adjacent-exon), stage7e (exon). Results under work/results/grid_grase
and work/results/grid_exon.
