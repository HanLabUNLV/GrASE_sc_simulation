#!/usr/bin/env bash
# stage7f: NB normalization test on the exon arm. norm=none already run (grid negbinom_EB);
# run offset_id ((1|id)+offset) and offset_only (offset, no RE -> no integration), format to
# the confirmed dexseq_file shape, score with the confirmed DEXSeq arm. Compare sens/spec.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"
COUNTS="$WORK/grid/counts/exon_sc.txt"
SF="$WORK/grid/counts/exon_sizefactors.txt"
OUT="$RESULTS_DIR/grid_exon"; mkdir -p "$OUT"
to_dexseq() { awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; print "groupID\tfeatureID\tpadj"; next}
  { ev=$h["event"]; g=$h["gene"]; bin=ev; sub(/^.*:/,"",bin); printf "%s\tE%s\t%s\n", g, bin, $h["padj"] }' "$1" > "$2"; }
for NM in offset_id offset_only; do
  tag="negbinom_${NM}"
  echo "[$(date '+%F %T')] $tag fit ..."
  "$RSCRIPT" "$GSIM/scripts/run_model.R" --counts "$COUNTS" --model negbinom --disp EB \
      --norm "$NM" --sizefactors "$SF" --out "$OUT/exon_${tag}.result.txt" --mc_cores 40 \
      > "$OUT/exon_${tag}.fit.log" 2>&1
  to_dexseq "$OUT/exon_${tag}.result.txt" "$OUT/exon_${tag}.dexseq.txt"
  mkdir -p "$OUT/eval_${tag}"
  "$RSCRIPT" "$GSIM/scripts/evaluate_tools_dexseq_vs_gt.R" \
      "$TRUTH_DIR/sim_exon_info" /dev/null /dev/null "$OUT/exon_${tag}.dexseq.txt" \
      "$OUT/eval_${tag}" "$DEDS_TRUTH" /nonexistent_jgt.txt > "$OUT/eval_${tag}.log" 2>&1 || true
  echo "[$(date '+%F %T')] $tag done"
done
echo "[$(date '+%F %T')] stage7f done"
