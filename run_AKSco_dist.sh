#!/bin/bash

#SBATCH -J Dist #Single job name for the entire JobArray

#SBATCH -o slurm/multi_%A_%a.out #standard output

#SBATCH -e slurm/multi_%A_%a.err #standard error

#SBATCH -p general #partition

#SBATCH -t 00:10:00 #running time

#SBATCH --mail-type=BEGIN

#SBATCH --mail-type=END

#SBATCH --mail-user=iancze@gmail.com

#SBATCH --mem-per-cpu 100 #memory request per node

#SBATCH -n 6

##SBATCH --cpus-per-task=10

#SBATCH --ntasks-per-node=2

python hostgen.py $SLURM_ARRAY_TASK_ID

julia --machinefile slurm/run${SLURM_ARRAY_TASK_ID}hosts.txt tests/parallel_test.jl
