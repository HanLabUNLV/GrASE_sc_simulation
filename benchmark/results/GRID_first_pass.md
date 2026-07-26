# scGrASE 3-counts x 3-models grid -- first pass (10M/80, 50bp)

jaquino's design realized on scARTist single-cell BAMs: three count representations
(SJ / exon / adjacent-exon) each run through the same three glmmTMB dispersion models
(beta-binomial / negative-binomial / mixed-binomial), so differences isolate the count
representation rather than each tool's native test. Upgrades over jaquino: (1) EB
dispersion moderation for all models via a single shared moderator; (2) the
mixed-binomial follows rMATS logit-normal logic (primary) with an EB variant; (3) a
single-fit + fixed-dispersion + Wald-contrast structure that removes jaquino's
asymmetric full/reduced LRT.

## Infrastructure (all committed)

GrASE_sc (package, branch singlecell):
- `R/model_negbinom.R` -- `test_model_negbinom` + `estimate_theta` (NB, EB theta).
- `R/model_mixedbinom.R` -- `test_model_mixedbinom` + `estimate_logitnormal_var` +
  `rmats_logitnormal_var` (mixed-binomial; disp = rmats primary / EB variant).
- Both reuse the shared `moderate_phi_log_scale` moderator and `run_one_comparison`
  scaffold unchanged; `exontest_functions.R` untouched (clean upstream merge).

GrASE_sc_simulation (application, branch singlecell):
- `scripts/run_model.R` -- thin model-agnostic driver (any model x any count type);
  comparison-aware EB; BH FDR + significance flag.
- `scripts/build_sj_counts.R` -- rMATS JCEC IJC/SJC -> sc (diff=inclusion, n=inc+skip).
- `scripts/build_exon_counts.R` -- DEXSeq bins -> sc (diff=bin, n=gene total); the
  scalable DEXSeq arm (no classic testForDEU).
- `scripts/build_adjexon_counts.R` -- GrASE bubble bipartitions -> sc (make_comparison,
  two comparisons/bubble).
- `scripts/evaluate_grid_sj.R` -- SJ arm vs junction GT (sens/spec/precision/recall).

Harness: `benchmark/stage7_grid.sh` builds the count tables and runs each cell.

## SJ arm: cross-check vs jaquino diff_params_4 (restricted universe, padj 0.05)

The SJ counts through the three models are the design-intended rMATS proxy.

| my SJ arm         | sens  | spec  |    | jaquino SJ (10M/80) | sens  | spec  |
|-------------------|-------|-------|----|---------------------|-------|-------|
| betabinom_EB      | 0.030 | 0.974 |    | betabinomial        | 0.293 | 0.930 |
| negbinom_EB       | 0.051 | 0.939 |    | negativebinomial    | 0.580 | 0.697 |
| mixedbinom_EB     | 0.152 | 0.848 |    | mixedbinomial       | 0.233 | 0.953 |
| mixedbinom_rmats  | 0.323 | 0.627 |    | (rMATS proxy)       |   -   |   -   |

Restricted universe = 4,192 grid-tested junction events (99 GT-positive). Full universe
= all 37,233 junction-GT events.

Internally coherent ordering: the rMATS-variance model (floored at 0.01) is the most
sensitive / least specific; EB moderation makes all models progressively more
conservative (mixed > nb > beta). But the numbers do NOT replicate jaquino cell-by-cell
-- notably his negative-binomial is his most sensitive SJ model (0.58) while ours is
near the bottom (0.05).

## Caveats / next (the model-accuracy + calibration audit)

The discrepancy is the expected subject of the accuracy audit -- candidate causes to
work through:

1. **Junction-GT universe mismatch.** `sim_junction_gt.txt` still spans 37,233 events
   (the full 94k-transcript annotation), while the filtered rMATS run (20,851 tx) tests
   4,193. The restricted universe controls for this, but the GT should be rebuilt on the
   filtered GTF for a clean full-universe number.
2. **EB over-shrinkage.** beta/nb sensitivity is near zero on SJ counts -- single-cell
   junction counts are highly overdispersed, so an EB-moderated dispersion suppresses
   power. jaquino used free per-feature dispersion; compare EB vs free on SJ.
3. **FDR + aggregation.** BH per contrast + min-padj-per-event may be harsher than
   jaquino's per-feature BH.
4. **Single-fit+contrast vs LRT.** the restructuring changes the power profile; verify
   against jaquino's helper output on a few features.
5. **Exon arm** (170k bins) not yet run -- the long pole; and the **adjacent-exon eval**
   (bubble -> exon-part GT) is still to be wired (reuse evaluate_bipartition_test.R).

Status: grid mechanics complete and validated end-to-end for SJ + adjacent-exon (8
cells); exon arm + adjexon/exon evaluation + the calibration audit are the follow-on.
