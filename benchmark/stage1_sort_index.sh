#!/usr/bin/env bash
# Stage 1: coordinate-sort + index scARTist genome BAMs (name-sorted as emitted).
# null -> group1 (control), deds -> group2 (case); rename to sample_NN.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"
mkdir -p "$BAMS_G1" "$BAMS_G2"

do_cell() {
  local arm=$1 n=$2 out
  case "$arm" in null) out=$BAMS_G1;; deds) out=$BAMS_G2;; esac
  local src="$READS_DIR/$arm/Bcell_genome.${n}.bam"
  local dst; dst="$out/sample_$(printf '%02d' "$n").bam"
  "$SAMTOOLS" sort -@ 2 -o "$dst" "$src"
  "$SAMTOOLS" index "$dst"
}
export -f do_cell
export SAMTOOLS BAMS_G1 BAMS_G2 READS_DIR

for arm in null deds; do
  for n in $(seq 1 "$NCELLS"); do echo "$arm $n"; done
done | xargs -P 24 -n 2 bash -c 'do_cell "$@"' _

echo "group1(null): $(ls "$BAMS_G1"/*.bam 2>/dev/null | wc -l)  group2(deds): $(ls "$BAMS_G2"/*.bam 2>/dev/null | wc -l)"
for d in "$BAMS_G1" "$BAMS_G2"; do
  for b in "$d"/*.bam; do "$SAMTOOLS" quickcheck "$b" || echo "QUICKCHECK FAIL: $b"; done
done
echo "quickcheck done"
