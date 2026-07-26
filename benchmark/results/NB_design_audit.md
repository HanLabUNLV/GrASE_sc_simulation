# Negative-binomial design audit: fixed `id` (jaquino) vs random `(1|id)` (mine)

Question: my adjacent-exon negbinom sensitivity (0.192) is well below jaquino's thesis
(0.346). Is it the random-`(1|id)`-vs-fixed-`id` design change?

## jaquino's exact model (adjexonnegbinomLRT.R)

  fit.full    <- glmmTMB(y ~ id + thot + groups:thot, family = nbinom2, se = FALSE)
  fit.reduced <- glmmTMB(y ~ id + thot,               family = nbinom2, se = FALSE)
  anova(fit.full, fit.reduced, test = "LRT")

`y = c(this, others)` stacked; `id` = cell as a FIXED factor; `thot` = this/others part;
`groups:thot` = the DTU interaction; test = LRT. Four differences from mine: (1) id fixed
vs random, (2) LRT vs Wald, (3) `id + thot` base vs `0 + groups:part`, (4) free theta vs
EB-fixed. This audit fixes theta free on both to isolate the id+test design.

## Head-to-head on the SAME adjacent-exon events (free theta both)

101 events where both models fit; betabinom-significant events used as high-confidence
truth.

|                              | jaquino (fixed id + LRT) | mine (random (1|id) + Wald) |
|------------------------------|--------------------------|-----------------------------|
| raw p<0.05 detected          | 30                       | 25                          |
| detects betabinom-sig (n=25) | 25 (100%)                | 23 (92%)                    |
| median p on those events     | 2.1e-17                  | 1.1e-8                      |

## Conclusion: a sensitivity/specificity tradeoff, driven by the cell effect

- Fixed `id` conditions out each cell's total EXACTLY (matched/conditional analysis),
  isolating the diff-vs-ref-by-group signal -> p-values ~9 orders of magnitude smaller on
  true events -> more survive FDR -> higher sensitivity (matches thesis 0.346).
- Random `(1|id)` PARTIALLY POOLS cell effects -> dilutes the signal -> less sensitive
  (0.192) but more specific (0.997 vs thesis 0.834). The extreme fixed-id p-values also
  produce more false positives (jaquino's lower spec).
- Fixed `id` is also LESS ROBUST: many events threw `rank-deficient conditional model`
  and failed (160 fixed nuisance params over 320 rows). That fragility is why the random
  effect was adopted.
- Conceptually, fixed-`id` NB ~ the beta-binomial (both condition on the cell total),
  which is why betabinom_EB already occupies the high-precision corner.

## Recommendation

The random-`(1|id)` + EB design is the deliberate high-specificity/robust choice; jaquino's
fixed-`id` + LRT is the high-sensitivity/less-robust one. Options:
- Keep random `(1|id)` (current) for specificity + robustness.
- Add jaquino's fixed-`id` + LRT as a selectable NB variant to reproduce the thesis
  sensitivity (accepting lower spec + convergence failures).
- Since fixed-`id` NB ~ beta-binomial, the betabinom_EB arm already provides the
  conditioned/high-sensitivity-per-fit behavior more stably; the NB arm's distinct value
  is the counts parameterization, best kept in its robust random-effect form.

Audit script: benchmark/analysis/nb_design_audit.R.
