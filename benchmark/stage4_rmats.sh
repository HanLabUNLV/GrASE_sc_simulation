#!/usr/bin/env bash
# Stage 4 (rMATS arm): rMATS-turbo on null(group1) vs deds(group2) sorted BAMs.
# One-shot (--task both), paired-end, readLength 50, annotation-only (no --novelSS,
# matching jaquino and the annotation-based v34 graph). Then map events to DEXSeq
# fragments with GrASE's map_rmats_graph.R.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"

POST=$RMATS_OUT/rmats_post_group1_group2
mkdir -p "$RMATS_OUT" "$POST" "$RMATS_OUT/tmp"

ls "$BAMS_G1"/sample_*.bam | paste -sd, > "$RMATS_OUT/b1.txt"   # null / control
ls "$BAMS_G2"/sample_*.bam | paste -sd, > "$RMATS_OUT/b2.txt"   # deds / case
echo "b1 (null): $(tr ',' '\n' < "$RMATS_OUT/b1.txt" | wc -l) bams | b2 (deds): $(tr ',' '\n' < "$RMATS_OUT/b2.txt" | wc -l) bams"

echo "[$(date '+%F %T')] rMATS prep+post..."
"$RMATS_PY" "$RMATS_SCRIPT" \
    --b1 "$RMATS_OUT/b1.txt" --b2 "$RMATS_OUT/b2.txt" \
    --gtf "$RMATS_GTF" -t paired --readLength "$READLEN" --nthread 40 \
    --od "$POST" --tmp "$RMATS_OUT/tmp" --task both

echo "[$(date '+%F %T')] rMATS done. JCEC files:"
ls "$POST"/*.MATS.JCEC.txt 2>/dev/null

echo "[$(date '+%F %T')] map rMATS events -> DEXSeq fragments (map_rmats_graph.R)..."
"$RSCRIPT" "$GRASE_PKG/Rpkg/scripts/map_rmats_graph.R" \
    "$GRAPHML_DIR" "$POST" "$RMATS_OUT/map_rmats"
echo "[$(date '+%F %T')] map done."
ls "$RMATS_OUT/map_rmats"/ 2>/dev/null | head
