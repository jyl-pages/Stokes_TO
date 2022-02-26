## Usage: python run_parallel_tasks.py

## Please set the following parameters
## The <COMMAND> will be executed for <REPEATED> times with maximum <MAX_PROC> processes
## For each process, it will be executed in a separate folder prefixed with <FOLDER_PREFIX>
COMMAND = "test_exe"  
REPEATED = 5  
MAX_PROC = 2  
FOLDER_PREFIX = "my_test_"  
## End setting parameters

# Note
# This script only works on Linux
# For Windows,  **(haven't been tested)**
  # change "os.wait" to "time.sleep(some_amount_of_time)" ;
  # change "if len(..." to "while len(..." .


import subprocess, os

command = "../" + COMMAND 
max_processes = MAX_PROC
repeated = REPEATED
prefix = FOLDER_PREFIX
processes = set()

for i in range(repeated):
    idx = i % max_processes
    folder = prefix + str(idx)
    if (i < max_processes):
        c = "mkdir " + folder + "; cd " + folder + "; " + command
    else :
        c = "cd " + folder + "; " + command
    print ("running " + c) 
    processes.add(subprocess.Popen(c, shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            [p for p in processes if p.poll() is not None])

for p in processes:
    if p.poll() is None:
        p.wait()


# Ref
# https://stackoverflow.com/questions/4992400/running-several-system-commands-in-parallel-in-python
# https://docs.python.org/3/library/subprocess.html