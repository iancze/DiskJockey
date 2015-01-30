#!/bin/bash

#SBATCH -J AKScoD #Single job name for the entire JobArray

#SBATCH -o slurm/AKScoD_%A_%a.out #standard output

#SBATCH -e slurm/AKScoD_%A_%a.err #standard error

#SBATCH -p general #partition

#SBATCH -t 01:00:00 #running time

#SBATCH --mail-type=BEGIN

#SBATCH --mail-type=END

#SBATCH --mail-user=iancze@gmail.com

#SBATCH --mem-per-cpu 500 #memory request per cpu

#SBATCH -n 51

##SBATCH --cpus-per-task=10
##SBATCH --ntasks-per-node=2

python hostgen.py $SLURM_ARRAY_TASK_ID

julia --machinefile slurm/run${SLURM_ARRAY_TASK_ID}hosts.txt burma_shave.jl -r $SLURM_ARRAY_TASK_ID scripts/AKSco.yaml --chain
