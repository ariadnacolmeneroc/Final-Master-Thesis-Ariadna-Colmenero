# -*- coding: utf-8 -*-
"""
From raw reads in FASTQ files to BAM files
Ariadna Colmenero Cobo de Guzman
Master's Final Thesis
2021-2022
Python Version: Python 3.9.7 
Code formating provider: Black
"""

"""
The following steps have been implemented using the Starlife cluster through a Linux Operating 
System (all jobs have been run using Slurm). First, using "ssh", we must remotely access the 
cluster through our user. Next, we access the folder where our files of interest are located. 
In this case, it is interesting to use "sshfs" to set up a remote file system to be able to 
see the results on our computer as we obtain them. It is very important to note that we need 
3 different files to launch each job. One must contain the python script with the command lines
(.py), another one contains the path to our files (.txt) and finally a .sh with the instructions
to execute the program on the cluster.
The follwoing lines of code have been executed for the step: Creation of an Interval List for the QC control metrics analysis.
"""


# First we need to import the different modules required:
import subprocess  # Allows us to spawn new processes, connect to their input/output/error pipes.
import argparse  # Helps us to parse those arguments that the program needs.
import time  # This module provides us with various time-related functions.
import sys  # It allows to operate on the interpreter as it provides access to variables and functions that interact strongly with the interpreter.

# For each of the selected files (indicated in the corresponding TXT with the paths), proceed as follows:
parser = argparse.ArgumentParser(
    prog="Interval_list.py", description="""From BED file to Interval list."""
)

parser.add_argument(
    "-o",
    "--outDir",
    dest="outDir",
    action="store",
    required=True,
    help="Path to out Working Directory.",
) # Path to the working directory where the BED file is located.

parser.add_argument(
    "-inp",
    "--inputfile",
    dest="inp",
    action="store",
    required=True,
    help="Path to the Bed file belonging to the analysed data.",
) # Arguments for the input BED file we want to transform.

parser.add_argument(
    "-sd",
    "--sd",
    dest="sd",
    action="store",
    required=True,
    help="Path to the sequence reference dictionary.",
) # The sequence dictionary, or BAM/VCF/IntervalList from which a dictionary can be extracted.


options = parser.parse_args()  # This will convert each argument to the appropriate type and then invoke the appropriate action.

# Define the different inputs as an object so that we do not have to enter their path each time they are used in the command line.
outDir = options.outDir
inp = options.inp
sd = options.sd

start_time = time.time() # At this point, the time it takes for the job to complete begins to count down.
 
LOG_FILE_OUT = open(outDir + "/log/log_out_il.txt", "w")  # In the error folder of the working directory, a file with the output log will be created.

LOG_FILE_ERR = open(outDir + "/log/log_err_il.txt", "w") # IMPORTANT CONSIDERATION: In the error folder of the working directory, a TXT is created that explains in which step the job is found. In case of error, it will indicate where it is in the code.

bashArguments = (
    "java -Xmx16G -jar /apps/PICARD/2.24.0/bin/picard.jar BedToIntervalList  I="
    + inp
    + " O=interval_list_TFM.interval_list SD="
    + sd
    + ""
) #  Through Picard's BedToIntervalList implementation (v2.24.0), we can obtain the .interval_lis output.

 # If there is an error in this pipeline, it shall be recorded in the error file indicated by "ERROR IN COMMAND:"
e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

if e != 0:
    LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
    LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
    sys.exit(1) # In case of an error in this sub-process, the job is stopped.
 
# If the sub-process has been completed, the log output file will show that the job has been completed and will also show how long it took.
 # If this is done, both error and log output files are closed and the process is finished.

LOG_FILE_OUT.write(
    "JOB DONE. Exectuion time in minutes: %s\n"
    % (round((time.time() - start_time) / 60, 4))
)
LOG_FILE_OUT.close()
LOG_FILE_ERR.close()
