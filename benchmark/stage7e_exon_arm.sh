#!/usr/bin/env bash
# stage7e: exon (DEXSeq-count) arm through the confirmed path. The exon bins are flat
# diff/n events, so the models run via run_model.R (confirmed modules: group_by_event,
# run_map/moderate, run_one_comparison + BH -- matching DEXSeq's own BH FDR). Each result
# is formatted into the confirmed dexseq_file shape (groupID/featureID=E+bin/padj) and
# scored by the CONFIRMED evaluate_tools_dexseq_vs_gt.R DEXSeq arm vs the same
# sim_exon_info exon-part GT (junction block skipped via a nonexistent junction_gt).
set -euo pipefail
source /mnt/data1/home/mirahan/scGrASE/benchmark/config.sh
COUNTS="$WORK/grid/counts/exon_sc.txt"
OUT="$RESULTS_DIR/grid_exon"; mkdir -p "$OUT"

to_dexseq() {  # <grid result> <dexseq file>
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; print "groupID\tfeatureID\tpadj"; next}
   { ev=$h["event"]; g=$h["gene"]; bin=ev; sub(/^.*:/,"",bin);
     printf "%s\tE%s\t%s\n", g, bin, $h["padj"] }' "$1" > "$2"
}

run_cell() {  # <model> <disp>
  local M="$1" D="$2" tag="${1}_${2}"
  echo "[$(date '+%F %T')] exon $tag: fit ..."
  "$RSCRIPT" "$GSIM/scripts/run_model.R" \
      --counts "$COUNTS" --model "$M" --disp "$D" \
      --out "$OUT/exon_${tag}.result.txt" --mc_cores 40 > "$OUT/exon_${tag}.fit.log" 2>&1
  to_dexseq "$OUT/exon_${tag}.result.txt" "$OUT/exon_${tag}.dexseq.txt"
  mkdir -p "$OUT/eval_${tag}"
  echo "[$(date '+%F %T')] exon $tag: eval (confirmed DEXSeq arm) ..."
  "$RSCRIPT" "$GSIM/scripts/evaluate_tools_dexseq_vs_gt.R" \
      "$TRUTH_DIR/sim_exon_info" /dev/null /dev/null "$OUT/exon_${tag}.dexseq.txt" \
      "$OUT/eval_${tag}" "$DEDS_TRUTH" /nonexistent_jgt.txt \
      > "$OUT/eval_${tag}.log" 2>&1 || echo "  (eval nonzero; junction arm skipped as expected)"
}

run_cell mixedbinom rmats
run_cell betabinom  EB
run_cell negbinom   EB
run_cell mixedbinom EB
echo "[$(date '+%F %T')] stage7e done -> $OUT"
