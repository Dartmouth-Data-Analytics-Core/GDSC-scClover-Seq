#!/bin/bash

#SBATCH --job-name=fbc_align
#SBATCH --nodes=1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.log
#========================================================#
# Example job: FBC1-only tRNA recovery (the pipeline's main mode).
# config.yaml's own defaults already match FBC1/fbc_only -- the --config
# flags below are restated explicitly anyway so this script never
# silently depends on config.yaml staying the same in the future.
#
# This job's own SBATCH header only needs to stay alive to coordinate --
# --profile cluster_profile submits each rule as its OWN separate sbatch
# job (see cluster_profile/config.yaml), sized to that rule's actual
# threads/maxtime/mem_mb (e.g. sc_tRNA_bowtie2 alone gets 16 cpus/60gb),
# instead of this parent job reserving one big allocation upfront sized
# for the worst-case combination of concurrently-running rules.
#========================================================#

set -euo pipefail

#----- Environment information
CONDA_BASE="/optnfs/common/miniconda3"
SNAKEMAKE_ENV="/dartfs/rc/nosnapshots/G/GMBSR_refs/envs/snakemake"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${SNAKEMAKE_ENV}"

# Scratch space for Snakemake's own per-rule conda envs and temp files --
# avoids filling up home-directory quota and known package-cache issues
# on this cluster (see CONTEXT.md).
export XDG_CACHE_HOME="/dartfs-hpc/scratch/${USER}/.cache"
export TMPDIR="/dartfs-hpc/scratch/${USER}/tmp"
export CONDA_PKGS_DIRS="/dartfs-hpc/scratch/${USER}/conda_pkgs"
mkdir -p "$XDG_CACHE_HOME" "$TMPDIR" "$CONDA_PKGS_DIRS"

# cluster_profile/config.yaml's own submission command writes each rule's
# SLURM log to slurm_logs/ -- must exist before the first job submits.
mkdir -p slurm_logs/

#----- LOGGER
cat <<EOF
#───────────────────────── Initialization ──────────────────────────#
Running Single-cell Clover-Seq (FBC1-only) with Snakemake $(snakemake --version)

Job:        $SLURM_JOB_NAME
Job ID:     $SLURM_JOB_ID
Node:       $(hostname)
Start time: $(date)
Work dir:   $(pwd)
Conda base: $CONDA_BASE
Snakemake:  $SNAKEMAKE_ENV
Binary:     $(which snakemake)
Profile:    cluster_profile/
#───────────────────────────────────────────────────────────────────#

SNAKEMAKE LOG:
EOF

#----- Invoke snakemake
# Always `snakemake -n` with these exact --config flags first, on your
# own copy of this script, before dropping -n and running it for real.
snakemake -s Snakefile \
	--conda-frontend conda \
	--use-conda \
	--profile cluster_profile \
	--rerun-incomplete \
	--keep-going \
	--config library_filter=FBC1 library_exclude=__none__ results_dir=fbc_only

#----- Capture exit status
SNAKEMAKE_EXIT=$?

#----- Final status
echo ""
echo "#------------------------ Job Complete ------------------------#"
echo "End time: $(date)"
if [ $SNAKEMAKE_EXIT -eq 0 ]; then
    echo "Status: SUCCESS"
else
    echo "Status: FAILED (exit code: $SNAKEMAKE_EXIT)"
fi

exit $SNAKEMAKE_EXIT
