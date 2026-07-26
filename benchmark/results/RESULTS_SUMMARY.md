# scGrASE single-cell DTU benchmark -- results summary

Config: scARTist single-cell simulated genome BAMs, 10M reads / 80 cells per arm
(null=group1, deds=group2), 50 bp reads. jaquino's 3-counts x 3-models DTU design,
realized on the confirmed GrASE pipeline with three upgrades: EB dispersion moderation
for all models (single shared moderator), the mixed-binomial follows rMATS logit-normal
logic, and a model-accuracy audit.

## Design and how replication is guaranteed

Three count representations x three glmmTMB dispersion models, so differences isolate the
count representation and the dispersion model, not each tool's native test:

| count type      | how counts are built          | how the model is run + scored                         |
|-----------------|-------------------------------|-------------------------------------------------------|
| adjacent-exon   | GrASE bubble bipartitions     | exontest.R (new --model branches) -> evaluate_bipartition_test.R |
| exon (DEXSeq)   | DEXSeq bin vs gene total      | run_model.R (flat, confirmed modules+BH) -> confirmed dexseq_file -> evaluate_tools DEXSeq arm |
| SJ (rMATS)      | rMATS IJC/SJC junction counts | run_model.R -> confirmed evaluate_tools junction arm  |

Models: betabinom (EB), negbinom (EB), mixedbinom (EB variant, and rmats variant =
rMATS logit-normal variance max(0.01,10*var(psi))). Every cell holds the confirmed
GT + evaluator + FDR constant and varies only the model, so numbers are comparable.

VALIDATION: the rebuilt grase reproduces the confirmed betabinom result exactly (identical
46-event significant set), and the adjacent-exon betabinom eval reproduces the confirmed
first-pass 0.785 precision / 0.216 recall (271 DTU genes) exactly. The comparison is legit.

## Main results -- ALL sim_type (nulls + DGE charged as FP), restricted universe, padj 0.05

sens = recall = TP/(TP+FN); spec = TN/(TN+FP); prec = TP/(TP+FP).

| count arm     | model            | sens  | spec  | prec  | TP   | FP    | FN   | TN     |
|---------------|------------------|-------|-------|-------|------|-------|------|--------|
| adjacent-exon | betabinom_EB     | 0.216 | 0.997 | 0.754 | 150  | 49    | 544  | 17429  |
| adjacent-exon | negbinom_EB      | 0.192 | 0.997 | 0.689 | 133  | 60    | 558  | 17392  |
| adjacent-exon | mixedbinom_EB    | 0.205 | 0.997 | 0.759 | 142  | 45    | 550  | 17432  |
| adjacent-exon | mixedbinom_rmats | 0.444 | 0.942 | 0.232 | 308  | 1022  | 386  | 16456  |
| exon          | betabinom_EB     | 0.415 | 0.991 | 0.574 | 1979 | 1471  | 2789 | 159863 |
| exon          | negbinom_EB      | 0.347 | 0.990 | 0.498 | 1653 | 1666  | 3115 | 159537 |
| exon          | mixedbinom_EB    | 0.357 | 0.992 | 0.560 | 1702 | 1340  | 3066 | 159990 |
| exon          | mixedbinom_rmats | 0.833 | 0.838 | 0.132 | 3974 | 26175 | 794  | 135158 |

DTU-only precision/recall (restricted, @0.05) for reference -- adjacent-exon: betabinom
0.785/0.216, negbinom 0.760/0.193, mixedbinom_EB 0.807/0.205, mixedbinom_rmats 0.627/0.444;
exon: betabinom 0.611/0.415, mixedbinom_rmats 0.481/0.834.

### Reading it
- Three EB models = high precision / modest recall; mixedbinom_rmats trades precision for
  recall and, under ALL (nulls charged), precision collapses (adj 0.23, exon 0.13 / 26k FP)
  -- driven by the tight 0.01 variance floor.
- Exon (DEXSeq-bin) counts are uniformly more sensitive than adjacent-exon (betabinom
  recall 0.42 vs 0.22) at lower precision -- consistent with jaquino's thesis.

## Comparison to jaquino thesis (10M/80, @0.05, sens/spec)

| arm      | model         | my sens | thesis sens | my spec | thesis spec |
|----------|---------------|---------|-------------|---------|-------------|
| adj-exon | betabinom     | 0.216   | 0.270       | 0.997   | 0.914       |
| adj-exon | mixedbinom_EB | 0.205   | 0.245       | 0.997   | 0.942       |
| adj-exon | negbinom      | 0.192   | 0.346       | 0.997   | 0.834       |
| exon     | betabinom     | 0.415   | 0.625       | 0.991   | 0.687       |

Comparable in regime; sensitivity same order, generally a bit lower. betabinom and
mixedbinom map closely. TWO caveats:
1. NEGBINOM diverges most (0.192 vs 0.346, adj) -- likely the audit item: my NB uses a
   random (1|id) cell effect vs jaquino's fixed id, which shifts power.
2. Specificity is NOT apples-to-apples: my restricted-universe ALL has a large null TN
   pool so spec sits ~0.99 vs the thesis 0.83-0.94 -- a universe/denominator difference,
   not a real model difference. (mixedbinom_rmats is a new rMATS-faithful variant not in
   the thesis.)

## rMATS native arm is degenerate on single-cell data

Native rMATS FDR is NOT a valid comparator here. Confirmed causal chain:
- native-rMATS recall 0.027 (13 TP) <- 99.0% of events get PValue == 1 <- rMATS sets p=1
  for any event with a ZERO-total-coverage sample (its zero-count filter) <- single-cell
  sparsity (80 cells/group, many zero-coverage per junction; 99.2% of events have >=1).
- Cross-tab: every SE event with a zero-coverage cell (2956/2956) -> PValue 1; all 11
  significant events have no zero cell. Example SE_1754: obvious dPSI 0.74, rMATS PValue=1.
- Controlled test on the matched junction GT: the gap is the TEST not the filters -- the
  proxy stays 23x more sensitive (recall 0.61 vs 0.027) even with rMATS's own filters on.
This is exactly why the design proxies rMATS with SJ-counts-through-glmmTMB. Detail:
RMATS_proxy_vs_native.md.

## Methods notes

- Environment: grase rebuilt into a local lib (work/env/Renviron.grase_sc prepends it,
  preserving jaquino:sylvia). check_grase_integrity.sh confirms it reproduces confirmed
  betabinom. Runners: system R 4.4.2 (jaquino 4.2 lib) for GrASE/eval; rmats.4.2.0 conda;
  smartSim python for HTSeq counting.
- NB fits are the compute bottleneck: estimate_theta fits count ~ 0+groups:part+(1|id)
  (observation-level random effect -> Laplace integration), on 2x-stacked rows, twice per
  bin (Poisson + nbinom2 for the identifiability LRT). Beta-binomial puts overdispersion in
  the family (closed form, no integration). See SESSION_SUMMARY.md for how to avoid it
  (moment/DESeq2-style theta; fixed cell effects; or note betabinom is the integration-free
  equivalent).

## Harness / provenance
- benchmark/: stage7b (adj-exon fit), stage7c (adj-exon eval), stage7d (SJ rMATS 3-way),
  stage7e (exon fit+eval), check_grase_integrity.sh.
- Committed to forks (branch singlecell, no AI signatures): GrASE_sc (model_negbinom.R,
  model_mixedbinom.R, exontest.R --model dispatch, NAMESPACE); GrASE_sc_simulation
  (run_model.R, build_{sj,exon,adjexon}_counts.R, evaluate_tools grid-padj override).
- Result tables: work/results/grid_grase/<model>/, work/results/grid_exon/,
  work/results/sj_rmats_compare/. Companion docs: GRID_all_metric.md, GRID_adjexon.md,
  RMATS_proxy_vs_native.md.

## Open items
- negbinom design audit: random (1|id) vs jaquino's fixed id (the 0.19->0.35 sens gap);
  and a cheaper moment/DESeq2-style theta to remove the integration cost.
- SJ arm sens/spec through the confirmed junction evaluator for the EB models (proxy done
  for mixedbinom_rmats).
- Plotting and more configs/depths (150 cells; 41M/83M).
