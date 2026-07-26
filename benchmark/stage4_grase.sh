#!/usr/bin/env bash
# Stage 4 (GrASE arm = adjacent-exon model, latest GrASE's OWN beta-binomial).
# Uses ~/DICE prebuilt filtered bipartitions (no split step). exoncnt aggregates
# per-cell DEXSeq counts onto the internal-AS splits; exontest runs the
# beta-binomial EB test (betabinom_EBapprox). null=group1, deds=group2.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"
mkdir -p "$GRASE_OUT"

# GrASE = internal alternative splicing + TSSTTS (alternative TSS/TTS). Run both;
# the -a flag is 'internal' or 'TSS' (the latter writes the TSSTTS files).
for ALT in internal TSSTTS; do
  aflag=$([ "$ALT" = internal ] && echo internal || echo TSS)
  echo "[$(date '+%F %T')] exoncnt ($ALT): aggregate counts onto bipartition splits..."
  "$RSCRIPT" "$GRASE_PKG/Rpkg/scripts/exoncnt.R" \
      -i "$BIPARTITION_DIR" \
      -o "$GRASE_OUT/bipartition.${ALT}.counts" \
      -c "$COUNTS_DIR" -t bipartition \
      --cond1 group1 --cond2 group2 -a "$aflag"

  echo "[$(date '+%F %T')] exontest ($ALT): beta-binomial EB..."
  "$RSCRIPT" "$GRASE_PKG/Rpkg/scripts/exontest.R" \
      --file "bipartition.${ALT}.exoncnt.combined.txt" \
      --outdir "$GRASE_OUT/bipartition.test" \
      --countdir "$GRASE_OUT/bipartition.${ALT}.counts/" \
      --splittype bipartition \
      --phi "phi.glmmtmb.${ALT}.txt" \
      --model betabinom_EBapprox \
      --cond1 group1 --cond2 group2 --mc_cores 40
done

echo "[$(date '+%F %T')] GrASE arm done. results:"
ls "$GRASE_OUT/bipartition.test"/test_bipartition.{internal,TSSTTS}_*.annotated.txt 2>/dev/null
