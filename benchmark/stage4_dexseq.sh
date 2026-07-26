#!/usr/bin/env bash
# Stage 4 (DEXSeq arm): classic DEXSeq DEU test on the per-cell exon counts.
# null=group1 (control) vs deds=group2 (case).
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"
mkdir -p "$WORK/DEXSeq/out"

echo "[$(date '+%F %T')] DEXSeq DEU test..."
"$RSCRIPT" "$GSIM/DEXSeq/dexseq.R" \
    --gff "$FLATTENED_GFF" \
    --cntdir "$COUNTS_DIR" \
    --outdir "$WORK/DEXSeq/out" \
    --cell1 group1 --cell2 group2

echo "[$(date '+%F %T')] DEXSeq arm done. outputs:"
ls "$WORK/DEXSeq/out"/ 2>/dev/null | head
