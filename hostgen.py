#!/usr/bin/env python

# Create a list of all the CPU's we've been allocated, in addition to the one
# currently running.
from subprocess import check_output
import os

check = lambda cmd: check_output(cmd, shell=True, universal_newlines=True)

slurm_job_nodelist = check("echo $SLURM_JOB_NODELIST")

hostlist = check("scontrol show hostname " + slurm_job_nodelist).split("\n")

hostname = check("hostname")[:-1]
print(hostname)

# hostlist = "holy2a03108\nholy2a03201\nholy2a03106".split("\n")

print(hostlist)

tasklist = check("echo $SLURM_JOB_CPUS_PER_NODE").split(",")
# tasklist = "4,2(x3),5".split(",")

# if two nodes have the same number of processes, they could be 5(x2),3,5(x3)
# expand this into a list that is the same length as hostlist
tasks = []
for task in tasklist:
    if "x" in task:
        ntask, ntimes = task.replace("(", "").replace(")", "").split("x")
        for i in range(int(ntimes)):
            tasks.append(int(ntask))
    else:
        tasks.append(int(task))

print(tasks)

hostfile = "hosts.txt"
if os.path.isfile(hostfile):
    os.remove(hostfile)

f = open(hostfile, "w")

for i, host in enumerate(hostlist):
    for j in range(tasks[i]):
        f.write(host + "\n")

f.close()

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
