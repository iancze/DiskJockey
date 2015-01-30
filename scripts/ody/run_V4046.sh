#!/bin/bash

#SBATCH -J V4046 #Single job name for the entire JobArray

#SBATCH -o slurm/V4046_%A_%a.out #standard output

#SBATCH -e slurm/V4046_%A_%a.err #standard error

#SBATCH -p general #partition

#SBATCH -t 12:00:00 #running time

#SBATCH --mail-type=BEGIN

#SBATCH --mail-type=END

#SBATCH --mail-user=iancze@gmail.com

#SBATCH --mem 16000 #memory request per node

#SBATCH -N 1 #ensure all jobs are on the same node

#SBATCH -n 26

julia burma_shave.jl -r $SLURM_ARRAY_TASK_ID scripts/V4046Sgr.yaml
