#!/usr/bin/env python

'''
The point of this script is to generate a `machinefile` for Julia. This is a
list of *additional* hostnames on which to start worker processes. That means
that if we want to have a total of 51 processes running (1 master and 50 workers)
then the hostfile will need to have 50 lines containing the hostnames of the 50
other nodes. This requires first figuring out all of the nodes which SLURM has allocated, how many cores (CPUs) have been allocated on each node and the hostname of the master process. Then we loop through and write out all of the hostnames except the current master host.
'''

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

hostlist = check("scontrol show hostname " + slurm_job_nodelist)[:-1].split("\n")

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

hostfile = "slurm/run{}hosts.txt".format(args.run)
if os.path.isfile(hostfile):
    os.remove(hostfile)

print("Hostlist", hostlist)
print("Tasks", tasks)

f = open(hostfile, "w")
skipped = False

for ntask, host in zip(tasks, hostlist):
    for j in range(ntask):
        if (host == masterhost) and (not skipped):
            skipped = True
        else:
            f.write(host + "\n")

f.close()
