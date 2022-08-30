# -*- coding: utf-8 -*-
"""
From raw reads in FASTQ files to BAM files
Ariadna Colmenero Cobo de Guzman
Master's Final Thesis
2021-2022
Python Version: Python 3.9.7 (default, Sep 16 2021, 16:59:28) [MSC v.1916 64 bit (AMD64)]
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

The follwoing lines of code have been executed for the step: From raw reads in FASTQ files to 
BAM files.

- INPUT FILES: FASTQ files containing the raw reads.
    
- OUTPUT FILES: 
"""

# First we need to import the different modules required:
import subprocess  # Allows us to spawn new processes, connect to their input/output/error pipes.
import argparse  # Helps us to parse those arguments that the program needs.
import time  # This module provides us with various time-related functions.
import sys  # It allows to operate on the interpreter as it provides access to variables and functions that interact strongly with the interpreter.

# For each of the selected files (indicated in the corresponding TXT with the paths), proceed as follows:

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="The following pipeline allows us to go from raw reads in FASTQ format to BAM files sorted and with the duplicates resulting from the previous PCR process removed."
    )  # Description of the programme shown below

    # Now,  we have to fill the ArgumentParser with information anout the program arguments:

    # First, we will need to define the working directory where, in this instance, the fastq files are located.
    parser.add_argument(  # Add a new argument to the job.
        "-o",  # Prefix used in the corresponding TXT to indicate the working directory.
        "--outDir",  # Alternative prefix used in the corresponding TXT to indicate the working directory.
        dest="outDir",  # Name of the attribute to be added to the object returned by the parse_args function.
        action="store",  # This indicates how the command-line arguments should be handled.
        required=True,  #  It is not optional. It is needed to run the job.
        help="Path to the working directory where the files in query are located.",  # Help option describing this argument.
    )

    # Next, we define the arguments for both FASTQ files (In the case of having only one, the argument referring to the second would be eliminated).

    parser.add_argument(
        "-fastq1",
        "--fastq1",
        dest="fastq1",
        action="store",
        help="Path to the folder where fastq1 is located in the directory.",
    )  # Add argument for the first FASTQ file.

    parser.add_argument(
        "-fastq2",
        "--fastq2",
        dest="fastq2",
        action="store",
        help="Path to the folder where fastq1 is located in the directory.",
    )  # Add argument for the second FASTQ file.

    parser.add_argument(
        "-refGenome",
        "--refGenome",
        dest="refGenome",
        action="store",
        help="Path to the folder where the reference genome is located.",
    )  # Add argument for reference genome. IMPORTANT CONSIDERATION: the genome must also be indexed in the same folder (.dict, .fai).

    options = (
        parser.parse_args()
    )  # This will convert each argument to the appropriate type and then invoke the appropriate action.
    # Define the different inputs as an object so that we do not have to enter their path each time they are used in the command line.
    outDir = options.outDir
    fastq1 = options.fastq1
    fastq2 = options.fastq2
    refGenome = options.refGenome

    start_time = (
        time.time()
    )  # At this point, the time it takes for the job to complete begins to count down.

    sample = fastq1.split("/")[1].split("_")[
        0
    ]  # Define the name of the samples (which is used for the files that are created as output) from the corresponding name in fastq1.

    LOG_FILE_OUT = open(
        outDir + "/error/log_out_" + sample + "_qc.txt", "w"
    )  # In the error folder of the working directory, a file with the output log will be created.

    LOG_FILE_ERR = open(
        outDir + "/error/log_err_" + sample + "_qc.txt", "w"
    )  # IMPORTANT CONSIDERATION: In the error folder of the working directory, a TXT is created that explains in which step the job is found. In case of error, it will indicate where it is in the code.

    bashArguments = (
        "bwa mem -t 8 -R '@RG\\tID:"
        + sample
        + "\\tLB:"
        + sample
        + "\\tSM:"
        + sample
        + "\\tPL:ILLUMINA' "
        + refGenome
        + " "
        + fastq1
        + "  > Bam/"
        + sample
        + ".sam"
    )  # Use of bwa mem to align raw reads against the reference genome. A SAM file with the alignment information is created as output with the same name.

    # If there is an error in this pipeline, it shall be recorded in the error file indicated by "ERROR IN COMMAND:"
    e = subprocess.call(
        bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR
    )

    if e != 0:
        LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
        LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
        sys.exit(1)  # In case of an error in this sub-process, the job is stopped.
    bashArguments = (
        "samtools view -b -u -h -@8 -o Bam/"  # Allows us to convert the SAM file to the corresponding BAM file.
        + sample
        + ".bam Bam/"
        + sample
        + ".sam"
    )  # Using the SAM file that has just been created as input, the corresponding BAM file is created using SAMTools.

    # The commands to create the error message for this particular thread are re-specified.
    e = subprocess.call(
        bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR
    )

    if e != 0:
        LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
        LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
        sys.exit(1)  # In case of an error in this sub-process, the job is stopped.
    bashArguments = (
        "samtools sort -@8 -o Bam_KIEL/"  # Se ordenan las alineaciones por las coordenadas más a la izquierda del archivo BAM. Así, se crea un output con el nombre *_sorted.bam.
        + sample
        + "_sorted.bam Bam_KIEL/"
        + sample
        + ".bam"
    )  # How to input the BAM created earlier is needed.

    # The error is then redefined for this sub-process.
    e = subprocess.call(
        bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR
    )

    if e != 0:
        LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
        LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
        sys.exit(1)  # In case of an error in this sub-process, the job is stopped.
    # This command lines index the coordinate-sorted BAM file For fast random access.
    bashArguments = "samtools index -@8 -b Bam_KIEL/" + sample + "_sorted.bam"

    # The error is then redefined for this sub-process.
    e = subprocess.call(
        bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR
    )

    if e != 0:
        LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
        LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
        sys.exit(1)  # In case of an error in this sub-process, the job is stopped.
    # The following thread, needs to indicate that it is done via Java and also the location where the Picard .jar is located.
    bashArguments = (
        "java -Xmx16G -jar /apps/PICARD/2.24.0/bin/picard.jar MarkDuplicates I= Bam/"  # Picard's MarkDuplicates utility is used.
        + sample
        + "_sorted.bam O=Bam/"  # Se necesita cómo introducir el archivo BAM ordenado e indexado.
        + sample
        + "_sorted_dup.bam M=Bam/"  # The output file to write marked records to.
        + sample
        + "_duplicates.txt"  # File to write duplication metrics to.
    )

    # The error is then redefined for this sub-process.
    e = subprocess.call(
        bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR
    )

    if e != 0:
        LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
        LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
        sys.exit(1)  # In case of an error in this sub-process, the job is stopped.
    # Use of SAMTools view as shown above.
    bashArguments = "samtools index -@8 -b Bam_KIEL/" + sample + "_sorted_dup.bam"

    # The error is then redefined for this sub-process.
    e = subprocess.call(
        bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR
    )

    if e != 0:
        LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
        LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
        sys.exit(1)  # In case of an error in this sub-process, the job is stopped.
    # Next, the CollectInsertSizeMetrics is implemented through Java useful metrics for validating library construction.
    bashArguments = (
        "java -Xmx16G -jar /apps/PICARD/2.24.0/bin/picard.jar CollectInsertSizeMetrics I=Bam_KIEL/"
        + sample
        + "_sorted_dup.bam O=Bam_KIEL/"  # O is the file to write the output to.
        + sample
        + "_insert_size_metrics.txt H=Bam_KIEL/"  # H is the file to write insert size Histogram chart to.
        + sample
        + "_insert_size_histogram.pdf M=0.5"  # Descart data with fewer than this percentage overall reads.
    )

    # The error is then redefined for this sub-process.
    e = subprocess.call(
        bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR
    )

    if e != 0:
        LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
        LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
        sys.exit(1)  # In case of an error in this sub-process, the job is stopped.
    # Finally, for reasons of space (and because they are not useful for subsequent processes), the following files are deleted:
    bashArguments = (
        "rm Bam_KIEL/"
        + sample
        + ".sam Bam_KIEL/"
        + sample
        + ".bam Bam_KIEL/"
        + sample
        + ".bam.bai"
    )

    # The error is then redefined for this sub-process.
    subprocess.call(bashArguments, shell=True)

    if e != 0:
        LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
        LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
        sys.exit(1)  # In case of an error in this sub-process, the job is stopped.
    # If all the sub-processes have been completed, the log output file will show that the job has been completed and will also show how long it took.
    # If this is done, both error and log output files are closed and the process is finished.
    LOG_FILE_OUT.write(
        "JOB DONE. Exectuion time in minutes: %s\n"
        % (round((time.time() - start_time) / 60, 4))
    )
    LOG_FILE_OUT.close()
    LOG_FILE_ERR.close()
