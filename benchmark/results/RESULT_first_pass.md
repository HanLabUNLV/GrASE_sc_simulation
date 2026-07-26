# scGrASE first-pass result (config: rd 10M, 80 cells, 50 bp reads)

DTU detection benchmark on scARTist single-cell simulated genome BAMs (null=group1
vs deds=group2, 80 cells each). Ground truth = scARTist injected DTU
(transcript_ds_status), rendered to exonic-part / junction GT via the (adapted)
GrASE_simulation truth builders. Evaluation = GrASE_simulation `evaluate_*` scripts.

## Detection arms
- GrASE (latest, own bubble beta-binomial `exontest`, model betabinom_EBapprox),
  internal + TSSTTS bipartitions, using ~/DICE prebuilt v34 partitions.
- rMATS turbo 4.2 (native FDR), paired, readLength 50, annotation-only.
- DEXSeq (classic testForDEU): DID NOT FINISH -- ran >4.8 h on 160 samples and was
  killed. Classic DEXSeq does not scale to single-cell sample counts. To be added
  via a scalable glmmTMB-based exon test (as jaquino's diff_params_4 did).

## Results (restricted universe = testable events; DTU genes; padj 0.05)

| method | precision | recall | F1 (micro) | macro F1 | TP | FP | FN |
|--------|-----------|--------|------------|----------|----|----|----|
| GrASE (internal+TSSTTS) | 0.63 | 0.34 | 0.45 | 0.66 | 1694 | 990 | 3240 |
| rMATS (native FDR)      | 0.93 | 0.028| 0.053| 0.91 | 14   | 1   | 495 |

GrASE internal-only @0.05: precision 0.785, recall 0.216 (271 testable DTU genes).

## Headline
On single-cell simulated data, GrASE's beta-binomial adjacent-exon/bubble model is
~12x more sensitive than rMATS's native test for DTU (recall 0.34 vs 0.028), at
lower precision (0.63 vs 0.93). rMATS is near-perfect precision but detects only
~15 of ~509 testable DTU events.

## Caveats / follow-ups
- rMATS here = its NATIVE FDR. jaquino's "rMATS SJ counts" arm instead re-tests the
  rMATS splice-junction counts with glmmTMB (more sensitive) -- a different method;
  add that arm for a like-for-like sensitivity comparison.
- DEXSeq arm pending a scalable test.
- One config (10M/80), read length 50, GrASE beta-binomial only. Other
  configs/depths/models are follow-on.
- Harness: ~/scGrASE/benchmark + modified clones on branch scartist-singlecell.
