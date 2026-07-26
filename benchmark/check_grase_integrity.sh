#!/usr/bin/env bash
# Integrity check: re-run betabinom_EBapprox on the SAME internal bipartition counts
# with the rebuilt (local-lib) grase, reusing the cached phi, and compare the
# significant calls + padj to the CONFIRMED output. If they match, the rebuilt grase
# reproduces jaquino's installed grase for betabinom -> safe to trust new-model numbers.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"
# use the rebuilt grase (local lib prepended, jaquino:sylvia preserved)
export R_ENVIRON_USER=/mnt/data1/home/mirahan/scGrASE/work/env/Renviron.grase_sc

CHECK=/mnt/data1/home/mirahan/scGrASE/work/grase_check
mkdir -p "$CHECK"
# reuse the cached phi so estimation is skipped (exercises moderation/run/adjust/annotate)
cp -f "$GRASE_OUT/bipartition.test/phi.glmmtmb.internal.txt" "$CHECK/phi.glmmtmb.internal.txt"

echo "[check] re-running betabinom_EBapprox with rebuilt grase (grase at:"
"$RSCRIPT" -e 'cat(find.package("grase")[1])'; echo ")"

"$RSCRIPT" "$GRASE_PKG/Rpkg/scripts/exontest.R" \
    --file "bipartition.internal.exoncnt.combined.txt" \
    --outdir "$CHECK" \
    --countdir "$GRASE_OUT/bipartition.internal.counts/" \
    --splittype bipartition \
    --phi "phi.glmmtmb.internal.txt" \
    --model betabinom_EBapprox \
    --cond1 group1 --cond2 group2 --mc_cores 40 > "$CHECK/run.log" 2>&1

NEW=$(ls "$CHECK"/test_bipartition.internal_betabinom_EBapprox.annotated.txt 2>/dev/null | head -1)
OLD="$GRASE_OUT/bipartition.test/test_bipartition.internal_betabinom_EBapprox.annotated.txt"
echo "[check] confirmed: $OLD"
echo "[check] rebuilt:   $NEW"
scount() { awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="significant")s=i} NR>1&&$s=="TRUE"{n++} END{print n}' "$1"; }
echo "[check] significant rows -- confirmed: $(scount "$OLD")  rebuilt: $(scount "$NEW")"
# compare the sorted (gene,event,contrast,padj,significant) projections
proj() { awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){h[$i]=i}} NR>1{printf "%s\t%s\t%s\t%.6g\t%s\n",$h["gene"],$h["event"],$h["contrast"],$h["padj"],$h["significant"]}' "$1" | sort; }
proj "$OLD" > "$CHECK/old.proj"; proj "$NEW" > "$CHECK/new.proj"
if diff -q "$CHECK/old.proj" "$CHECK/new.proj" >/dev/null; then
  echo "[check] PASS: padj + significance identical across all events"
else
  echo "[check] DIFF: $(diff "$CHECK/old.proj" "$CHECK/new.proj" | grep -c '^<') confirmed-only, $(diff "$CHECK/old.proj" "$CHECK/new.proj" | grep -c '^>') rebuilt-only rows"
  diff "$CHECK/old.proj" "$CHECK/new.proj" | head -12
fi
