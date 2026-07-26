#!/usr/bin/env bash
# stage7b: run the NEW adjacent-exon models (negbinom_EB, mixedbinom_EB,
# mixedbinom_rmats) through the SAME exontest.R + same internal bipartition counts +
# same annotation as the confirmed beta-binomial arm. Only --model changes. Output goes
# to the same bipartition.test dir (distinct per-model filenames), so
# evaluate_bipartition_test.R scores them exactly like betabinom vs the same sim_exon_info GT.
set -euo pipefail
source /mnt/data1/home/mirahan/scGrASE/benchmark/config.sh
export R_ENVIRON_USER=/mnt/data1/home/mirahan/scGrASE/work/env/Renviron.grase_sc   # rebuilt grase

TESTDIR="$GRASE_OUT/bipartition.test"
for MODEL in negbinom_EB mixedbinom_EB mixedbinom_rmats; do
  echo "[$(date '+%F %T')] exontest ($MODEL) on internal counts..."
  "$RSCRIPT" "$GRASE_PKG/Rpkg/scripts/exontest.R" \
      --file "bipartition.internal.exoncnt.combined.txt" \
      --outdir "$TESTDIR" \
      --countdir "$GRASE_OUT/bipartition.internal.counts/" \
      --splittype bipartition \
      --phi "phi.${MODEL}.internal.txt" \
      --model "$MODEL" \
      --cond1 group1 --cond2 group2 --mc_cores 40 \
      > "$TESTDIR/run_${MODEL}.log" 2>&1
  echo "[$(date '+%F %T')]   -> $(ls "$TESTDIR"/test_bipartition.internal_${MODEL}.annotated.txt 2>/dev/null || echo FAILED)"
done
echo "[$(date '+%F %T')] stage7b done."
