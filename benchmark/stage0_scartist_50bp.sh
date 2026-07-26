#!/usr/bin/env bash
# Stage 0: regenerate scARTist reads at 50bp for BOTH arms (null+deds), config 10M/80.
# Reuses cached groundtruth + sim_counts (read_len-independent); rebuilds pbsim3(deds)
# and all reads at 50bp with Emp50 profiles via config.benchmark50.yaml.
set -euo pipefail
cd /mnt/data1/home/mirahan/scRNAsim/scARTist
SMK=/mnt/data1/home/mirahan/miniconda3/envs/beers2/bin/snakemake
PROJ_ROOT=/mnt/data1/home/mirahan/scRNAsim
CORES=80

rm -rf results/reads/rd_10000000_cells_80          # drop stale 125bp null reads
TARGETS=""
for arm in null deds; do
  for i in $(seq 1 80); do
    TARGETS="$TARGETS results/reads/rd_10000000_cells_80/${arm}/Bcell_genome.${i}.bam"
  done
done

echo "[$(date '+%F %T')] Stage0: 160 BAMs @50bp (Emp50), ${CORES} cores"
# NOTE: --configfile takes nargs='+' in snakemake 7.x, so it must come FIRST
# (immediately followed by an option) or it swallows the positional BAM targets.
"${SMK}" --configfile config/config.benchmark50.yaml \
    --cores "${CORES}" --use-singularity \
    --singularity-args "--bind ${PROJ_ROOT}" \
    ${TARGETS}
echo "[$(date '+%F %T')] done"
echo "null BAMs: $(ls results/reads/rd_10000000_cells_80/null/*.bam 2>/dev/null | wc -l)"
echo "deds BAMs: $(ls results/reads/rd_10000000_cells_80/deds/*.bam 2>/dev/null | wc -l)"
