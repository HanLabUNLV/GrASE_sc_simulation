# Adjacent-exon grid (3 models) -- CONFIRMED-pipeline results (10M/80, 50bp, internal)

The three grid models run through the SAME exontest.R driver + same internal bipartition
counts + same FDR (nested_BH + independent filtering) + same annotation as the confirmed
beta-binomial arm, and scored with the SAME evaluate_bipartition_test.R against the SAME
sim_exon_info exon-part GT. Only --model varies, so the numbers are directly comparable.

Environment integrity: the rebuilt (local-lib) grase reproduces the confirmed betabinom
internal result EXACTLY -- identical 46-event significant set, padj to ~5-6 sig figs.
And betabinom's eval below (0.785 / 0.216, 271 genes) reproduces the confirmed first-pass
number exactly -> the comparison is valid.

## Restricted DTU, padj 0.05 (micro precision/recall)

| model              | precision | recall | TP  | FP  | FN  |
|--------------------|-----------|--------|-----|-----|-----|
| betabinom_EBapprox | 0.785     | 0.216  | 150 | 41  | 544 |
| negbinom_EB        | 0.760     | 0.193  | 133 | 42  | 558 |
| mixedbinom_EB      | 0.807     | 0.205  | 142 | 34  | 550 |
| mixedbinom_rmats   | 0.627     | 0.444  | 308 | 183 | 386 |

271 testable DTU genes. Event-level significant counts (padj<0.01): betabinom 46,
negbinom 47, mixedbinom_EB 40, mixedbinom_rmats 216.

The three EB models cluster (prec 0.76-0.81, recall 0.19-0.22), consistent with jaquino's
adjacent-exon regime (sens 0.24-0.35). mixedbinom_rmats is 2x more sensitive / less
precise -- the liberal floored-variance behavior.

## Open: mixed-rMATS proxy (liberal) vs native rMATS (conservative)

The confirmed end-to-end NATIVE rMATS arm was ultra-conservative (recall 0.028). The
mixedbinom_rmats proxy is liberal (recall 0.444 here). Two confounds are stacked and must
be separated:
1. SUBSTRATE: native rMATS scored JUNCTION counts vs junction GT; the 0.444 above is
   mixedbinom_rmats on ADJACENT-EXON (bubble) counts vs exon-part GT. Different data.
2. TEST: native rMATS's own likelihood + FDR + filters (avg count >= 10, PSI in
   [0.05,0.95], |dPSI| cutoff) vs the proxy (binomial GLMM, logit-normal variance fixed
   at max(0.01,10*var(psi)), Wald + BH). The 0.01 variance floor is very tight given 160
   cells -> high power.

To isolate the TEST effect: run mixedbinom_rmats on the SJ counts and score it with the
CONFIRMED junction evaluator (evaluate_tools_dexseq_vs_gt.R) + the MATCHED junction GT
(sim_junction_gt_matched.txt, 5,279), then compare to native rMATS FDR on the same
junctions. That is the still-to-wire SJ arm (grid padj through the confirmed evaluator).

## Method (all confirmed code)

- exontest.R --model {negbinom_EB,mixedbinom_EB,mixedbinom_rmats} (new dispatch branches;
  reuse run_map_for_comparison / moderate_phi_log_scale / run_one_comparison /
  adjust_pvalues / add_significant / annotation).
- benchmark/stage7b_grid_adjexon.sh (run) + stage7c_eval_adjexon.sh (eval).
- grase rebuilt into local lib (work/env/Renviron.grase_sc prepends it, preserving
  jaquino:sylvia); check_grase_integrity.sh validates it reproduces confirmed betabinom.
