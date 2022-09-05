# -*- coding: utf-8 -*-
"""
Variant Calling - Mutect2
Ariadna Colmenero Cobo de Guzman
Master's Final Thesis
2021-2022
"""

########################## Variant Calling - Mutect2 ##########################################

"""
The following steps have been implemented using the Starlife cluster through a Linux Operating 
System (all jobs have been run using Slurm). First, using "ssh", we must remotely access the 
cluster through our user. Next, we access the folder where our files of interest are located. 
In this case, it is interesting to use "sshfs" to set up a remote file system to be able to 
see the results on our computer as we obtain them. It is very important to note that we need 
3 different files to launch each job. One must contain the python script with the command lines
(.py), another one contains the path to our files (.txt) and finally a .sh with the instructions
to execute the program on the cluster. Moreover, we have used Black as the code formatter. 

The follwoing lines of code were executed for the step: Application of the Variant Calling 
Mutect2 to the sorted and marked duplicated BAM files (Obtained through the FASTQTOBAM.PY).

IMPORTANT CONSIDERATION: This pipeline can be used for both targeted and whole exome analysis.
                         However, it is necessary to change the genomic intervals. These, are
                         indicated in the code.
"""

# First we need to import the different modules required:
import subprocess  # Allows us to spawn new processes, connect to their input/output/error pipes.
import argparse  # Helps us to parse those arguments that the program needs.
import time  # This module provides us with various time-related functions.
import sys  # It allows to operate on the interpreter as it provides access to variables and functions that interact strongly with the interpreter.

# First, we define what our programme is and add a brief description of it.
parser = argparse.ArgumentParser(
    prog="VARIANT_CALLER_MUTECT2.py",
    description="""Variant calling code implementation.""", #
)

# For each of the selected files (indicated in the corresponding TXT with the paths), proceed as follows:
parser.add_argument(  # Add a new argument to the job.
        "-o",  # Prefix used in the corresponding TXT to indicate the working directory.
        "--outDir",  # Alternative prefix used in the corresponding TXT to indicate the working directory.
        dest="outDir",  # Name of the attribute to be added to the object returned by the parse_args function.
        action="store",  # This indicates how the command-line arguments should be handled.
        required=True,  #  It is not optional. It is needed to run the job.
        help="Path to the working directory where the files in query are located.",  # Help option describing this argument.
    )

# Next, we define the arguments for the BAM file as we have just performed with the directory.
parser.add_argument(
    "-bam",
    "--bam",
    dest="bam",
    action="store",
    help="bam file generated in the Fastq_Bam_AC2.py",
)

# Then, we define the differnt intputs so that we do not have to enter their path each time they are used in the command line.
options = parser.parse_args()
outDir = options.outDir
bam = options.bam

start_time = time.time() # At this point, the time it takes for the job to complete begins to count.

sample = bam.split("/")[1].replace(".bam", "").split("_")[0] # Define the name of the samples (which is used for the files that are created as output) from the corresponding name in the BAM file.

LOG_FILE_OUT = open(outDir + "/log_VC/log_out_" + sample + ".txt", "w") # In the error folder of the working directory, a file with the output log will be created.

LOG_FILE_ERR = open(outDir + "/log_VC/log_err_" + sample + ".txt", "w") # IMPORTANT CONSIDERATION: In the error folder of the working directory, a TXT is created that explains in which step the job is found. In case of error, it will indicate where it is in the code.

######
### BASERECALIBRAOR
######

# In the next step, we must change the genomic intervals depending on whether our samples are Target Sequencing or Whole Exome Sequencing. 
# In this case, for Target sequencing samples we use the BED file belonging to the panel shown in Supplementary Figure 1 (Master's Thesis report) with the 168 genes related to lymphomagenesis. 
# On the other hand, for the 5 WES samples that passed the quality control, the Exome V6 is used.
bashArguments = (
    "java -Xmx16G -jar /apps/GATK/4.1.8.1/gatk-package-4.1.8.1-local.jar BaseRecalibrator -I " # Execute BaseRecalibrator from GATK's .jar.
    + bam # Refers to the BAM file as the input file.
    + " -L Info_VC/BEDfile_Mutations_Panel_Design_no11q_168gens.bed -O Mutect2/" # If we work with WES, we must include in the intervals argument the corresponding exome version here.
    + sample # At this point the first output is generated: a TXT with information about the recalibration.
    + "_recal.txt -R Ref_Genome/genome.fa --known-sites Info_VC/dbsnp_138.hg19.vcf --known-sites Info_VC/1000G_phase1.indels.hg19.vcf" # Both files for known sites can be used for both sequencing approaches.
) # Implementation of BaseRecalibrator (GATK 4.1.8.1).

e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
# If there is an error in this pipeline, it shall be recorded in the error file indicated by "ERROR IN COMMAND:"
if e != 0:
    LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
    LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
    sys.exit(1)# In case of an error in this sub-process, the job is stopped.

######
### APPLYBQSR
######

bashArguments = (
    "java -Xmx16G -jar /apps/GATK/4.1.8.1/gatk-package-4.1.8.1-local.jar ApplyBQSR -I " # Acess to the ApllyBQSR tool.
    + bam # We indicate which is the input file (the same as in the previous step).
    + " -L Info_VC/BEDfile_Mutations_Panel_Design_no11q_168gens.bed -O Mutect2/" # Here, we have to insert the corresponding BED file for the genomic intervals.
    + sample # A BAM file containing the recalibrated read data.
    + "_BQSR.bam -R Ref_Genome/genome.fa --bqsr-recal-file Mutect2/"
    + sample
    + "_recal.txt" # Here, it needs the TXT file with the recalibration information created in the previous step via BaseRecalibrator. 
) # Implementation of ApplyBQSR (GATK 4.1.8.1)

# If there is an error in this pipeline, it shall be recorded in the error file indicated by "ERROR IN COMMAND:"
if e != 0:
    LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
    LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
    sys.exit(1) # In case of an error in this sub-process, the job is stopped.

######
### MUTECT2
######

# Here, we can call somatic short mutations via MUTECT2. For this, we have to proceed:
bashArguments = (
    "java -Xmx16G -jar /apps/GATK/4.1.8.1/gatk-package-4.1.8.1-local.jar Mutect2 -I Mutect2/" # Path to the Mutect2 tool from GATK 4.1.8.1.
    + sample # Use the BAM file created via APLYBQSR  as the input file.
    + "_BQSR.bam -L Info_VC/BEDfile_Mutations_Panel_Design_no11q_168gens.bed -O Mutect2/" # Reference to the genomic intervals (change it if necessary).
    + sample
    + ".vcf -R Ref_Genome/genome.fa -tumor " # Creation of the VCF file which is going to be the output.
    + sample
    + " --native-pair-hmm-threads 8 --max-reads-per-alignment-start 0" # How many threads should a native pairHMM implementation use.
                                                                       # Maximum number of reads to retain per alignment start position. 
)

# If there is an error in this pipeline, it shall be recorded in the error file indicated by "ERROR IN COMMAND:"
e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

if e != 0:
    LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
    LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
    sys.exit(1) # In case of an error in this sub-process, the job is stopped.
    
######
### FILTERMUTECTCALLS
######

# Finally, we can filter variants in a Mutect2 VCF callset:
bashArguments = (
    "java -Xmx16G -jar /apps/GATK/4.1.8.1/gatk-package-4.1.8.1-local.jar FilterMutectCalls -V Mutect2/" # We indicate the path to the FilterMutectCalls tool.
    + sample
    + ".vcf -O Mutect2/" # We use as the input file, the VCF obtained through the Mutect2 tool.
    + sample
    + "_filt.vcf --native-pair-hmm-threads 8 --reference Ref_Genome/genome.fa"
) # A TXT is retrieved with the filters applied to the raw output for Mutect2.

e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

# If there is an error in this pipeline, it shall be recorded in the error file indicated by "ERROR IN COMMAND:"
if e != 0:
    LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
    LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
    sys.exit(1) # In case of an error in this sub-process, the job is stopped.

 # If all the sub-processes have been completed, the log output file will show that the job has been completed and will also show how long it took.
    # If this is done, both error and log output files are closed and the process is finished.
    LOG_FILE_OUT.write(
        "JOB DONE. Exectuion time in minutes: %s\n"
        % (round((time.time() - start_time) / 60, 4))
    )
    LOG_FILE_OUT.close()
    LOG_FILE_ERR.close()
