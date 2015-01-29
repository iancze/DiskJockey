#!/bin/bash

hostlist=$(scontrol show hostname $SLURM_JOB_NODELIST)
rm -f hosts.txt

for f in $hostlist
do
  echo $f >> hosts.txt
done

echo $SLURM_JOB_CPUS_PER_NODE >> hosts.txt
echo $SLURM_TASKS_PER_NODE >> hosts.txt
