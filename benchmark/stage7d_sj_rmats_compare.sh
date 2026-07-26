#!/usr/bin/env bash
# stage7d: isolate the mixed-rMATS-proxy vs native-rMATS test difference on JUNCTIONS.
# Same CONFIRMED junction evaluator, same MATCHED junction GT (5,279), same universe.
# Three runs:
#   native       -- native rMATS FDR (baseline conservative)
#   proxy_filt   -- mixedbinom_rmats padj, rMATS count/PSI filters ON  (isolates the test)
#   proxy_nofilt -- mixedbinom_rmats padj, filters OFF                 (isolates the filters)
set -euo pipefail
source /mnt/data1/home/mirahan/scGrASE/benchmark/config.sh

GT_EXON="$TRUTH_DIR/sim_exon_info"
JGT="$TRUTH_DIR/sim_junction_gt_matched.txt"       # correct, filtered-run-matched GT
RMATS_POST="$RMATS_OUT/rmats_post_group1_group2"
MAPDIR="$RMATS_OUT/map_rmats"
DUMMY="/dev/null"                                   # dexseq arm not used here
PROXY="$WORK/grid/results/sj__mixedbinom_rmats.txt"
OUT="$RESULTS_DIR/sj_rmats_compare"; mkdir -p "$OUT"

run() {  # <subdir> <extra args...>
  local sub="$1"; shift
  mkdir -p "$OUT/$sub"
  "$RSCRIPT" "$GSIM/scripts/evaluate_tools_dexseq_vs_gt.R" \
      "$GT_EXON" "$RMATS_POST" "$MAPDIR" "$DUMMY" "$OUT/$sub" "$DEDS_TRUTH" "$JGT" "$@" \
      > "$OUT/${sub}.log" 2>&1 || echo "[warn] $sub returned nonzero (see $OUT/${sub}.log)"
}

echo "[$(date '+%F %T')] native rMATS (matched GT)..."
run native
echo "[$(date '+%F %T')] proxy + rMATS filters ON..."
run proxy_filt   "$PROXY" mixedbinom_rmats yes
echo "[$(date '+%F %T')] proxy + filters OFF..."
run proxy_nofilt "$PROXY" mixedbinom_rmats no
echo "[$(date '+%F %T')] stage7d done -> $OUT"
