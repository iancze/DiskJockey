#!/bin/bash

hostlist=$(scontrol show hostname $SLURM_JOB_NODELIST)
rm -f hosts

for f in $hostlist
do
  echo $f >> hosts
done
