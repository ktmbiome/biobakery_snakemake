#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=192g
#SBATCH --cpus-per-task 48
#SBATCH --gres=lscratch:800

source myconda
mamba activate biobakery_basic

export CONFIG_FILE="config.yaml"

snakemake --cores $SLURM_CPUS_PER_TASK \
    --max-jobs-per-second 1 --max-status-checks-per-second 0.01 \
    --latency-wait 120 --configfile $CONFIG_FILE
