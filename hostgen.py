#!/usr/bin/env python

# Figure out what job index we submitted with, if any
import argparse
parser = argparse.ArgumentParser(description="Create a file of hostnames.")
parser.add_argument("run", type=int, default=0, help="JobArray index.")
args = parser.parse_args()

# Create a list of all the CPU's we've been allocated, in addition to the one
# currently running.
from subprocess import check_output
import os

check = lambda cmd: check_output(cmd, shell=True, universal_newlines=True)

slurm_job_nodelist = check("echo $SLURM_JOB_NODELIST")

hostlist = check("scontrol show hostname " + slurm_job_nodelist).split("\n")

masterhost = check("hostname").split(".")[0]

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

hostfile = "slurm/run{}hosts.txt".format(args["run"])
if os.path.isfile(hostfile):
    os.remove(hostfile)

f = open(hostfile, "w")
skipped = False

for i, host in enumerate(hostlist):
    for j in range(tasks[i]):
        if (host == masterhost) and (not skipped):
            skipped = True
        else:
            f.write(host + "\n")

f.close()
