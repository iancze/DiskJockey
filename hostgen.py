#!/usr/bin/env python

# Create a list of all the CPU's we've been allocated, in addition to the one
# currently running.
from subprocess import check_output
import os

slurm_job_nodelist = check_output("echo $SLURM_JOB_NODELIST", shell=True, universal_newlines=True)

print(slurm_job_nodelist)

hostlist = check_output("scontrol show hostname " + slurm_job_nodelist, shell=True, universal_newlines=True)

print(hostlist)

hostfile = "hosts.txt"
os.remove(hostfile)

# hostlist=$(scontrol show hostname $SLURM_JOB_NODELIST)
# rm -f hosts.txt
#
# # The way SLURM/Julia works, we might need to try to figure out our current
# # host, and make sure we don't add that to the list of hosts to SSH to.
#
# for f in $hostlist
# do
#   for i in {1..5}
#   do
#     echo $f >> hosts.txt
#   done
# done

#echo $SLURM_JOB_CPUS_PER_NODE >> hosts.txt
#echo $SLURM_TASKS_PER_NODE >> hosts.txt
