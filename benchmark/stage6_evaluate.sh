#!/usr/bin/env bash
# Stage 6: evaluate all three arms against the scARTist-derived ground truth.
#   - GrASE (bipartition beta-binomial) vs exonic-part GT  -> evaluate_bipartition_test.R
#   - rMATS (junction GT) + DEXSeq (exon-bin GT)           -> evaluate_tools_dexseq_vs_gt.R
# TP/FP/TN/FN + precision/recall/F1 at padj {0.01,0.05,0.1,0.2}, full + restricted.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"
mkdir -p "$RESULTS_DIR/grase" "$RESULTS_DIR/tools"

# --- GrASE arm ---
GRASE_TEST=$(ls "$GRASE_OUT"/bipartition.test/test_bipartition.internal_*.annotated.txt 2>/dev/null | paste -sd,)
echo "[$(date '+%F %T')] GrASE eval; test files: $GRASE_TEST"
"$RSCRIPT" "$GSIM/scripts/evaluate_bipartition_test.R" \
    "$GRASE_TEST" \
    "$TRUTH_DIR/sim_exon_info" \
    "$RESULTS_DIR/grase" \
    "$DEDS_TRUTH"

# --- rMATS + DEXSeq arms ---
DEXSEQ_FILE=$(ls "$WORK/DEXSeq/out"/*.txt 2>/dev/null | head -1)
echo "[$(date '+%F %T')] rMATS+DEXSeq eval; dexseq file: $DEXSEQ_FILE"
"$RSCRIPT" "$GSIM/scripts/evaluate_tools_dexseq_vs_gt.R" \
    "$TRUTH_DIR/sim_exon_info" \
    "$RMATS_OUT/rmats_post_group1_group2" \
    "$RMATS_OUT/map_rmats" \
    "$DEXSEQ_FILE" \
    "$RESULTS_DIR/tools" \
    "$DEDS_TRUTH" \
    "$TRUTH_DIR/sim_junction_gt.txt"

echo "[$(date '+%F %T')] evaluation done. summary tables:"
ls "$RESULTS_DIR"/grase/*summary* "$RESULTS_DIR"/tools/*summary* 2>/dev/null
