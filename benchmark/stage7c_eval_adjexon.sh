#!/usr/bin/env bash
# stage7c: evaluate the adjacent-exon grid (betabinom baseline + 3 new models) with the
# CONFIRMED evaluate_bipartition_test.R against the SAME sim_exon_info exon-part GT.
# Every model used the same counts, FDR, and annotation -> only the model differs, so
# the sens/precision are directly comparable (and betabinom must match its confirmed run).
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"

TESTDIR="$GRASE_OUT/bipartition.test"
OUT="$RESULTS_DIR/grid_grase"
mkdir -p "$OUT"

for MODEL in betabinom_EBapprox negbinom_EB mixedbinom_EB mixedbinom_rmats; do
  TF="$TESTDIR/test_bipartition.internal_${MODEL}.annotated.txt"
  [ -f "$TF" ] || { echo "[eval] MISSING $TF"; continue; }
  mkdir -p "$OUT/$MODEL"
  echo "[$(date '+%F %T')] eval $MODEL ..."
  "$RSCRIPT" "$GSIM/scripts/evaluate_bipartition_test.R" \
      "$TF" "$TRUTH_DIR/sim_exon_info" "$OUT/$MODEL" "$DEDS_TRUTH" \
      > "$OUT/${MODEL}.eval.log" 2>&1 || echo "[eval] $MODEL FAILED (see $OUT/${MODEL}.eval.log)"
done
echo "[$(date '+%F %T')] stage7c done. summaries in $OUT/<model>/"
