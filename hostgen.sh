#!/bin/bash

hostlist=$(scontrol show hostname $SLURM_JOB_NODELIST)
rm -f hosts.txt

# The way SLURM/Julia works, we might need to try to figure out our current
# host, and make sure we don't add that to the list of hosts to SSH to.

for f in $hostlist
do
  for i in {1..5}
  do
    echo $f >> hosts.txt
  done
done

#echo $SLURM_JOB_CPUS_PER_NODE >> hosts.txt
#echo $SLURM_TASKS_PER_NODE >> hosts.txt
