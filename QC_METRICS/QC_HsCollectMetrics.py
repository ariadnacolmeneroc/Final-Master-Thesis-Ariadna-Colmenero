# -*- coding: utf-8 -*-
"""
From raw reads in FASTQ files to BAM files
Ariadna Colmenero Cobo de Guzman
Master's Final Thesis
2021-2022
Python Version: Python 3.9.7 
Code formating provider: Black
"""

########################## From Fastq file to Bam files pipeline  #############################

"""
The following steps have been implemented using the Starlife cluster through a Linux Operating 
System (all jobs have been run using Slurm). First, using "ssh", we must remotely access the 
cluster through our user. Next, we access the folder where our files of interest are located. 
In this case, it is interesting to use "sshfs" to set up a remote file system to be able to 
see the results on our computer as we obtain them. It is very important to note that we need 
3 different files to launch each job. One must contain the python script with the command lines
(.py), another one contains the path to our files (.txt) and finally a .sh with the instructions
to execute the program on the cluster.
The follwoing lines of code have been executed for the step: Quality Control for BAM files. 
"""

# First we need to import the different modules required:
import subprocess  # Allows us to spawn new processes, connect to their input/output/error pipes.
import argparse  # Helps us to parse those arguments that the program needs.
import time  # This module provides us with various time-related functions.
import sys  # It allows to operate on the interpreter as it provides access to variables and functions that interact strongly with the interpreter.

# For each of the selected files (indicated in the corresponding TXT with the paths), proceed as follows:
parser = argparse.ArgumentParser(
    prog="QC_HsCollect.py", description="""BAM FILES QUALITY CONTROL"""
)

# First, we will need to define the working directory where, in this instance, the fastq files are located.
parser.add_argument(  # Add a new argument to the job.
    "-o",  # Prefix used in the corresponding TXT to indicate the working directory.
    "--outDir",  # Alternative prefix used in the corresponding TXT to indicate the working directory.
    dest="outDir",  # Name of the attribute to be added to the object returned by the parse_args function.
    action="store",  # This indicates how the command-line arguments should be handled.
    required=True,  #  It is not optional. It is needed to run the job.
    help="Path to the working directory where the files in query are located.",  # Help option describing this argument.
    )

parser.add_argument(
    "-bam",
    "--bamfile",
    dest="bam",
    action="store",
    required=True,
    help="Path to the BAM for which we want to test the Quality.",
) # Argument for the BAM file.

parser.add_argument(
    "-ref",
    "--refGenome",
    dest="ref",
    action="store",
    required=True,
    help="Path to the reference genome in a fasta format.",
)

parser.add_argument(
    "-bait",
    "--bait_intervals",
    dest="bait",
    action="store",
    required=True,
    help="Path to the Bait Intervals in a interval_list format.",
) # Path to the BAIT which is going to be the same as the TARGET file.

# This will convert each argument to the appropriate type and then invoke the appropriate action.
options = parser.parse_args()

outDir = options.outDir
bam = options.bam
bait = options.bait

start_time = time.time() # At this point, the time it takes for the job to complete begins to count down.

sample = bam.split("/")[-1].split("_")[0] # Define the name of the samples (which is used for the files that are created as output).

LOG_FILE_OUT = open(outDir + "/log/log_out_" + sample + ".txt", "w")  # In the error folder of the working directory, a file with the output log will be created.

LOG_FILE_ERR = open(outDir + "/log/log_err_" + sample + ".txt", "w") # IMPORTANT CONSIDERATION: In the error folder of the working directory, a TXT is created that explains in which step the job is found. In case of error, it will indicate where it is in the code.

bashArguments = (
    "java -Xmx16G -jar /apps/PICARD/2.24.0/bin/picard.jar CollectHsMetrics  I="
    + bam
    + " O=output_hs_metrics_"
    + sample
    + ".txt R="
    + ref
    + " BAIT_INTERVALS="
    + bait
    + " TARGET_INTERVALS="
    + bait
    + ""
) # Through the CollectHsMetrics implementation we can obtain a TXT with the metrics.

# If there is an error in this pipeline, it shall be recorded in the error file indicated by "ERROR IN COMMAND:"
e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

if e != 0:
    LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
    LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
    sys.exit(1) # In case of an error in this sub-process, the job is stopped.

# If all the sub-process has been completed, the log output file will show that the job has been completed and will also show how long it took.
    # If this is done, both error and log output files are closed and the process is finished.
LOG_FILE_OUT.write(
        "JOB DONE. Exectuion time in minutes: %s\n"
        % (round((time.time() - start_time) / 60, 4))
    )
LOG_FILE_OUT.close()
LOG_FILE_ERR.close()

