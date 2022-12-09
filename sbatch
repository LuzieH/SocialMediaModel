#!/bin/sh
#SBATCH --partition=big
#SBATCH --cpus-per-task=48
#SBATCH --time=7-0 
#SBATCH --job-name="socpde"

export OPENBLAS_NUM_THREADS=$SLURM_CPUS_ON_NODE
julia --project=. -L dev.jl -e "runslurm(a)"
