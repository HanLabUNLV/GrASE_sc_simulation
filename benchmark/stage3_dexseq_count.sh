#!/usr/bin/env bash
# Stage 3: per-cell DEXSeq exon counting with the v34 bygene flattened GFF.
# -s no (unstranded) matches jaquino's diff_params_4 counting. Output cleaned of
# HTSeq '_' summary lines. Count files -> count_files/{group1,group2}/sample_NN_counts.txt.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"
mkdir -p "$COUNTS_DIR/group1" "$COUNTS_DIR/group2"

count_one() {
  local grp=$1 bam=$2 base out
  base=$(basename "$bam" .bam)
  out="$COUNTS_DIR/$grp/${base}_counts.txt"
  "$HTSEQ_PY" "$DEXSEQ_COUNT_PY" \
      --format bam --order pos --paired yes --stranded no \
      "$FLATTENED_GFF" "$bam" "$out"
  grep -v '^_' "$out" | sed 's/"//g' > "$out.tmp" && mv "$out.tmp" "$out"
}
export -f count_one
export HTSEQ_PY DEXSEQ_COUNT_PY FLATTENED_GFF COUNTS_DIR

for b in "$BAMS_G1"/*.bam; do echo "group1 $b"; done > /tmp/scg_count_list.$$
for b in "$BAMS_G2"/*.bam; do echo "group2 $b"; done >> /tmp/scg_count_list.$$
xargs -P 24 -a /tmp/scg_count_list.$$ -n 2 bash -c 'count_one "$@"' _
rm -f /tmp/scg_count_list.$$

echo "group1 counts: $(ls "$COUNTS_DIR"/group1/*_counts.txt 2>/dev/null | wc -l)  group2: $(ls "$COUNTS_DIR"/group2/*_counts.txt 2>/dev/null | wc -l)"
echo "sample row (group1 sample_01):"; head -2 "$COUNTS_DIR"/group1/sample_01_counts.txt 2>/dev/null
