#!/bin/bash

#SBATCH --job-name=OCMarray
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-2
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --partition=standard
#SBATCH --mail-type=END,FAIL
#========================================================#
# Example job: cellranger multi, one array task per well.
#
# This is this pipeline's own PREREQUISITE, not part of the Snakemake
# pipeline itself -- run this (or your own equivalent) first if
# `cellranger multi` hasn't already been run for your wells. The
# Snakemake pipeline never invokes cellranger; it only ever reads this
# run's own CSV and outs/ directory (config.yaml's multi_csv/multi_out --
# see README's "Configuration" section).
#
# One array task per well: task 1 reads OCM1.csv and writes to
# 236860_1/ (--id), task 2 reads OCM2.csv and writes to 236860_2/. Run
# from the directory holding OCM1.csv/OCM2.csv (see the README's "Well
# configuration" example for what one of those CSVs contains). Extend
# --array and add matching OCM<N>.csv files for more wells; add your own
# --id/--csv naming as needed. Add --account=... below if your
# allocation requires one.
#========================================================#

set -euo pipefail

module load cellranger

cellranger multi \
    --id=236860_${SLURM_ARRAY_TASK_ID} \
    --csv=OCM${SLURM_ARRAY_TASK_ID}.csv \
    --localmem=128 \
    --localcores=16
