#!/usr/bin/env bash
# scGrASE benchmark configuration: paths, interpreters, and the single-cell config.
# Sourced by the stage scripts. Heterogeneous runners (see plan): GrASE/DEXSeq-R
# in system R, rMATS in its conda env, dexseq_count.py in an HTSeq env.

# ---- config (first pass) ---------------------------------------------------
export RD=10000000
export NCELLS=80
export READLEN=50           # matches jaquino + scARTist reads regenerated at 50bp

# ---- interpreters ----------------------------------------------------------
# System R 4.4.2 with the full Bioc + glmmTMB + grase stack (jaquino's 4.2 lib).
# Set R_LIBS explicitly so the harness is HOME-independent (does not rely on ~/.Renviron).
export RSCRIPT=/usr/bin/Rscript
export R_LIBS=/mnt/data1/home/jaquino/R/x86_64-pc-linux-gnu-library/4.2

# rMATS turbo 4.2.0
export RMATS_PY=/mnt/data1/home/mirahan/miniconda3/envs/rmats.4.2.0/bin/python
export RMATS_SCRIPT=/mnt/data1/home/mirahan/miniconda3/envs/rmats.4.2.0/rmats_turbo_v4_2_0/rmats.py

# dexseq_count.py needs HTSeq+pysam (smartSim env); pair with system samtools.
export HTSEQ_PY=/mnt/data1/home/mirahan/miniconda3/envs/smartSim/bin/python
export DEXSEQ_COUNT_PY=/mnt/data1/home/jaquino/R/x86_64-pc-linux-gnu-library/4.2/DEXSeq/python_scripts/dexseq_count.py
export SAMTOOLS=/usr/bin/samtools

# ---- reference data (all mutually consistent v34, from ~/DICE) -------------
export GRAPHML_DIR=/mnt/data1/home/mirahan/DICE/graphml.v34            # per-gene .graphml
export PERGENE_GFF_DIR=/mnt/data1/home/mirahan/DICE/dexseq.gff         # per-gene dexseq.gff
export BIPARTITION_DIR=/mnt/data1/home/mirahan/DICE/bipartition.filtered
export FLATTENED_GFF=/mnt/storage/jaquino/scRNAseq_sim_pt2/grase/graphml.dexseq.v34/gencode.v34.dexseq.bygene.gff
# GTF filtered to scARTist's 20,851 simulated transcripts, so rMATS's fromGTF event
# universe matches the simulation (avoids the flood of spurious empty-side events
# from the full 94k-transcript annotation). Built by stage: see work/ref/.
export RMATS_GTF=/mnt/data1/home/mirahan/scGrASE/work/ref/gencode.v34.scartist20851.gtf

# ---- scARTist outputs (the substrate) --------------------------------------
export SCARTIST=/mnt/data1/home/mirahan/scRNAsim/scARTist
export SCARTIST_ART_SIF=${SCARTIST}/resources/art.sif
export READS_DIR=${SCARTIST}/results/reads/rd_${RD}_cells_${NCELLS}   # null/ and deds/ BAMs
export DEDS_TRUTH=${SCARTIST}/results/sim_counts/rd_${RD}_cells_${NCELLS}/simulation_deds.txt
export NULL_TRUTH=${SCARTIST}/results/sim_counts/rd_${RD}_cells_${NCELLS}/simulation_null.txt

# ---- harness clones + work dirs --------------------------------------------
export SCGRASE=/mnt/data1/home/mirahan/scGrASE
export GSIM=${SCGRASE}/GrASE_simulation          # cloned+modified benchmark scripts
export GRASE_PKG=${SCGRASE}/GrASE                # cloned GrASE (Rpkg scripts)
export WORK=${SCGRASE}/work
export BAMS_G1=${WORK}/bams/group1               # null (control)
export BAMS_G2=${WORK}/bams/group2               # deds (case)
export COUNTS_DIR=${WORK}/DEXSeq/count_files     # group1/ group2/
export RMATS_OUT=${WORK}/rMATS
export GRASE_OUT=${WORK}/grase
export TRUTH_DIR=${WORK}/truth
export RESULTS_DIR=${WORK}/results
