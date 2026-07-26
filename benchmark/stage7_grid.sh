#!/usr/bin/env bash
# stage7_grid.sh -- the 3-counts x 3-models DTU grid (jaquino's design, on scARTist).
# Builds the three count tables (sj / exon / adjexon) in the shared sc format, then runs
# each (count type) x (model + dispersion) cell through the package models via the
# application driver run_model.R. Server-specific paths come from config.sh; the count
# builders and the driver are the committed GrASE_sc_simulation scripts.
#
# Grid cells per count type:
#   betabinom  disp EB      (beta-binomial, EB-moderated)
#   negbinom   disp EB      (negative-binomial, EB-moderated)
#   mixedbinom disp rmats   (rMATS proxy, primary -- rMATS logit-normal variance)
#   mixedbinom disp EB      (mixed-binomial, EB-moderated variance -- variant)
#
# Usage: bash stage7_grid.sh [counts]   counts = sj,exon,adjexon (default: sj,adjexon;
#   exon is the long pole -- 170k bins -- run it explicitly / in the background).

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "${HERE}/config.sh"

SCRIPTS="${GSIM}/scripts"
GRID_COUNTS="${WORK}/grid/counts"
GRID_RES="${WORK}/grid/results"
mkdir -p "${GRID_COUNTS}" "${GRID_RES}"

WHICH="${1:-sj,adjexon}"
MC="${MC_CORES:-32}"
export R_LIBS

run_r() { "${RSCRIPT}" "$@"; }

## ---- 1. build the count tables ----------------------------------------------
build_sj() {
  local out="${GRID_COUNTS}/sj_sc.txt"
  [ -s "$out" ] && { echo "[grid] sj counts exist"; return; }
  echo "[grid] building sj counts..."
  run_r "${SCRIPTS}/build_sj_counts.R" \
    --rmats_dir "${RMATS_OUT}/rmats_post_group1_group2" --out "$out"
}
build_exon() {
  local out="${GRID_COUNTS}/exon_sc.txt"
  [ -s "$out" ] && { echo "[grid] exon counts exist"; return; }
  echo "[grid] building exon counts..."
  run_r "${SCRIPTS}/build_exon_counts.R" \
    --group1_dir "${COUNTS_DIR}/group1" --group2_dir "${COUNTS_DIR}/group2" --out "$out"
}
build_adjexon() {
  local out="${GRID_COUNTS}/adjexon_sc.txt"
  [ -s "$out" ] && { echo "[grid] adjexon counts exist"; return; }
  echo "[grid] building adjexon counts..."
  run_r "${SCRIPTS}/build_adjexon_counts.R" \
    --counts "${GRASE_OUT}/bipartition.internal.counts/bipartition.internal.exoncnt.combined.txt" \
    --out "$out"
}

## ---- 2. run one grid cell ---------------------------------------------------
run_cell() {
  local count_type="$1" model="$2" disp="$3"
  local counts="${GRID_COUNTS}/${count_type}_sc.txt"
  local out="${GRID_RES}/${count_type}__${model}_${disp}.txt"
  [ -s "$counts" ] || { echo "[grid] MISSING counts: $counts"; return 1; }
  echo "[grid] === ${count_type} x ${model}/${disp} ==="
  run_r "${SCRIPTS}/run_model.R" \
    --counts "$counts" --model "$model" --disp "$disp" \
    --out "$out" --mc_cores "$MC" --padj_threshold 0.05
}

## ---- 3. orchestrate ---------------------------------------------------------
IFS=',' read -ra TYPES <<< "$WHICH"
for ct in "${TYPES[@]}"; do
  case "$ct" in
    sj)      build_sj ;;
    exon)    build_exon ;;
    adjexon) build_adjexon ;;
    *) echo "[grid] unknown count type: $ct"; exit 1 ;;
  esac
done

for ct in "${TYPES[@]}"; do
  run_cell "$ct" betabinom  EB
  run_cell "$ct" negbinom   EB
  run_cell "$ct" mixedbinom rmats
  run_cell "$ct" mixedbinom EB
done

echo "[grid] done. results -> ${GRID_RES}"
ls -la "${GRID_RES}"
