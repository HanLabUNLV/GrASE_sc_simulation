# scGrASE benchmark runbook -- full command sequence

The harness that orchestrates the single-cell DTU benchmark (scARTist BAMs -> GrASE /
rMATS / DEXSeq count arms -> 3-counts x 3-models glmmTMB grid -> evaluation). Each stage
sources `config.sh`. Results/findings for each phase are in `results/*.md`.

IMPORTANT: `config.sh` holds SERVER-SPECIFIC absolute paths and interpreter locations
(system R + jaquino's R_LIBS, rmats.4.2.0 conda, smartSim python, ~/DICE reference,
scARTist output dirs). Edit it for any other environment before running. Stages run from
the `~/scGrASE` workspace (paths are relative to it).

## 0. Config
    # edit paths/interpreters first
    benchmark/config.sh              # sourced by every stage; not run directly

## 1. Base pipeline (substrate + native arms + GT + first-pass eval)
    bash benchmark/stage0_scartist_50bp.sh    # scARTist: regen null+deds BAMs at 50bp (10M/80)
    bash benchmark/stage1_sort_index.sh       # coord-sort + index -> group1=null, group2=deds
    bash benchmark/stage3_dexseq_count.sh     # per-cell DEXSeq bin counts (160 cells)
    bash benchmark/stage4_rmats.sh            # rMATS-turbo on the FILTERED 20,851-tx GTF
    bash benchmark/stage4_grase.sh            # GrASE bubble beta-binomial (exoncnt.R + exontest.R)
    bash benchmark/stage4_dexseq.sh           # classic DEXSeq testForDEU (deferred -- does not scale)
    # ground truth (GrASE_simulation scripts, scARTist-adapted):
    #   scripts/infer_diff_exons_gt.R -> work/truth/sim_exon_info      (exon-part GT)
    #   scripts/infer_junctions_gt.R  -> work/truth/sim_junction_gt.txt + sim_junction_gt_matched.txt
    bash benchmark/stage6_evaluate.sh         # first-pass eval: GrASE + rMATS(native) + DEXSeq

## 2. Rebuild grase with the new models + integrity check
    # NAMESPACE already exports the new fns; install a superset into a LOCAL lib so it
    # shadows jaquino's install non-destructively:
    R_LIBS=~/scGrASE/Rlib:<jaquino lib> \
      R CMD INSTALL --no-multiarch --library=~/scGrASE/Rlib GrASE/Rpkg
    # a custom Renviron prepends the local lib, preserving jaquino:sylvia:
    #   work/env/Renviron.grase_sc
    bash benchmark/check_grase_integrity.sh   # confirm rebuilt grase reproduces confirmed betabinom

## 3. Grid -- adjacent-exon arm (through the CONFIRMED exontest.R pipeline)
    bash benchmark/stage7b_grid_adjexon.sh    # exontest.R --model {negbinom_EB,mixedbinom_EB,mixedbinom_rmats}
    bash benchmark/stage7c_eval_adjexon.sh    # evaluate_bipartition_test.R vs sim_exon_info

## 4. Grid -- exon (DEXSeq-count) arm
    # exon counts (flat bins) are built + run via run_model.R, formatted to the confirmed
    # dexseq_file shape, scored by the confirmed evaluate_tools DEXSeq arm:
    bash benchmark/stage7e_exon_arm.sh        # build_exon_counts.R -> run_model.R (3 models) -> eval

## 5. Grid -- SJ (rMATS-proxy) arm + native-vs-proxy investigation
    # SJ counts: scripts/build_sj_counts.R -> work/grid/counts/sj_sc.txt
    #   Rscript scripts/run_model.R --counts sj_sc.txt --model mixedbinom --disp rmats \
    #     --out work/grid/results/sj__mixedbinom_rmats.txt
    bash benchmark/stage7d_sj_rmats_compare.sh  # confirmed junction evaluator 3x (native / proxy+filter / proxy-nofilter)
    Rscript benchmark/analysis/rmats_discrepancy_trace.R   # why native rMATS = PValue 1 (zero-count filter)

## 6. NB normalization test (size-factor offset)
    Rscript benchmark/analysis/compute_exon_sizefactors.R  # -> work/grid/counts/exon_sizefactors.txt
    bash benchmark/stage7f_nb_norm.sh          # negbinom norm={offset_id,offset_only} vs baseline

## 7. NB design audit (fixed id vs random (1|id))
    Rscript benchmark/analysis/nb_design_audit.R           # jaquino fixed-id+LRT vs random-id+Wald

## Superseded / provenance only
- `stage7_grid.sh` -- the FIRST grid pass (SJ + adjacent-exon via run_model.R + plain BH +
  a wrong-GT evaluator). Replaced for adjacent-exon by stage7b/7c (exontest.R + confirmed
  eval). Kept for provenance; do not use for the adjacent-exon numbers.

## Result docs (results/)
RESULTS_SUMMARY.md (top-level), GRID_all_metric.md, GRID_adjexon.md, RMATS_proxy_vs_native.md,
NB_design_audit.md, NB_normalization_test.md, plus the first-pass RESULT_first_pass.md /
GRID_first_pass.md.
