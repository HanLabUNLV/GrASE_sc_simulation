# Why the mixed-rMATS proxy is liberal vs native rMATS conservative

Question: the confirmed native rMATS arm is ultra-conservative (recall ~0.027) while the
mixedbinom_rmats proxy is liberal. How much of the gap is the TEST vs the FILTERS vs the
SUBSTRATE?

## Controlled comparison (all confounds removed except the one under test)

Same JUNCTION counts, SAME matched junction GT (sim_junction_gt_matched.txt, 5,279),
SAME confirmed evaluator (evaluate_tools_dexseq_vs_gt.R), SAME restricted universe. The
only change is the significance source (native FDR vs proxy padj) and, for the last row,
whether rMATS's own count/PSI filters are applied to the proxy.

RESTRICTED, DTU sim_type, padj 0.05:

| run                         | TP  | FP  | FN  | precision | recall |
|-----------------------------|-----|-----|-----|-----------|--------|
| native rMATS FDR            |  13 |   1 | 470 | 0.929     | 0.027  |
| proxy + rMATS filters ON    | 296 |  71 | 187 | 0.807     | 0.613  |
| proxy, filters OFF          | 331 | 117 | 152 | 0.739     | 0.685  |

## Conclusion: it's the TEST, not the filters

- With rMATS's OWN filters applied, the proxy still jumps from recall 0.027 -> 0.613
  (a 23x increase). So the conservatism of native rMATS is NOT its count/PSI filtering.
- Removing the filters adds only ~0.07 recall (0.613 -> 0.685) at some precision cost.
  Filters are a minor contributor.
- Therefore native rMATS's conservatism lives in its actual test/likelihood.

## Where the test differs (audit targets)

The proxy = binomial GLMM, logit-normal variance FIXED at max(0.01, 10*var(psi)), Wald
test of the group contrast, BH. Native rMATS = paired constrained-vs-unconstrained
likelihood-ratio test with its own variance handling and a |dPSI| > c null cutoff. Prime
suspects for the gap, to check next:
1. NULL hypothesis: rMATS tests H0 |dPSI| <= c (composite null with a cutoff), the proxy
   tests H0 dPSI = 0 (point null). A cutoff null is much more conservative.
2. Variance magnitude actually used: compare the proxy's max(0.01, 10*var(psi)) per event
   to rMATS's fitted/regularized variance on the same events -- if rMATS's is larger, that
   alone drives conservatism.
3. Wald vs LRT test statistic, and FDR over rMATS's tested set vs the grid's set.

## Tracking down a discrepancy example -> native rMATS is DEGENERATE here

Traced one event where the proxy calls DTU but native rMATS does not (SE_1754,
ENSG00000115053.16, GT-positive DTU, GT dPSI 0.817):

  counts:  g1 incl 30.3 / skip 112  (PSI 0.209)   g2 incl 325 / skip 12  (PSI 0.951)
  observed dPSI = 0.742  (huge, obvious)
  native rMATS: IncLevelDifference = -0.780  BUT PValue = 1, FDR = 1  -> not significant
  proxy:        z^2 = 2011, padj = 0, var used = 0.293  -> significant (correct)

An obvious dPSI of 0.74 with PValue EXACTLY 1 is contradictory. Checking the whole run:

  native rMATS JCEC events: 4193
  PValue == exactly 1.0 : 4150  (99.0%)
  PValue < 0.05         : 33
  events with |IncLevelDifference| >= 0.3 : 356
    of those PValue == 1 : 354 (99.4%);  PValue < 0.05 : 2

=> Native rMATS returns PValue = 1 for ~99% of events, including nearly every obvious
large-change event. This is NOT principled conservatism -- the rMATS statistical test is
DEGENERATE on this single-cell data (p defaults to 1 almost everywhere). The confirmed
native-rMATS arm's recall 0.027 reflects that degeneracy, not a conservative model. The
proxy is detecting real DTU that native rMATS is simply blind to.

ROOT CAUSE (confirmed): rMATS-turbo sets PValue = 1 for any event with a ZERO-total-
coverage sample (its zero-count filter). Cross-tab over SE events (PValue==1 vs "has a
zero-coverage cell"):

              has_zero_cell
  PValue==1   FALSE  TRUE
    FALSE        19     0
    TRUE          5  2956

EVERY event with a zero-coverage cell (2956/2956) -> PValue 1; all 11 significant events
have NO zero cell. 99.2% of SE events have >=1 zero-coverage cell (SE_1754: 17/160). So:

  native-rMATS recall 0.027  <-  99% PValue==1  <-  zero-count filter  <-  single-cell
  sparsity (80 cells/group, many zero-coverage per junction).

Native rMATS is structurally incompatible with single-cell sparsity: its zero-count
filter assumes bulk replicates that all have coverage. This is exactly WHY jaquino's
design proxies rMATS with SJ-counts-through-glmmTMB rather than native rMATS -- the mixed
model handles zero-coverage cells and recovers the real DTU (e.g. SE_1754) that native
rMATS discards. Takeaway: native rMATS FDR is not a valid comparator on single-cell data;
the glmmTMB SJ proxy is the design-correct way to bring rMATS's junction representation
to single cell.

## How it was run
benchmark/stage7d_sj_rmats_compare.sh -> the confirmed evaluator 3x (native, proxy_filt,
proxy_nofilt) via the new optional grid-padj override (backward compatible; native path
unchanged). SJ proxy padj = work/grid/results/sj__mixedbinom_rmats.txt.
