# -*- coding: utf-8 -*-
"""
Variant Calling - Mutect2
Ariadna Colmenero Cobo de Guzman
Master's Final Thesis
2021-2022
"""

########################## Variant Calling - Mutect2 #############################

"""
The following steps have been implemented using the Starlife cluster through a Linux Operating 
System (all jobs have been run using Slurm). First, using "ssh", we must remotely access the 
cluster through our user. Next, we access the folder where our files of interest are located. 
In this case, it is interesting to use "sshfs" to set up a remote file system to be able to 
see the results on our computer as we obtain them. It is very important to note that we need 
3 different files to launch each job. One must contain the python script with the command lines
(.py), another one contains the path to our files (.txt) and finally a .sh with the instructions
to execute the program on the cluster.

The follwoing lines of code have been executed for the step: Application of the Variant Calling 
Mutect2 to the sorted and marked duplicated BAM files. 

IMPORTANT CONSIDERATION: This pipeline can be used for both targeted and whole exome analysis.
                         However, it is necessary to change some of the files. All the indications
                         are listed in the following code.
"""

# First we need to import the different modules required:
import subprocess  # Allows us to spawn new processes, connect to their input/output/error pipes.
import argparse  # Helps us to parse those arguments that the program needs.
import time  # This module provides us with various time-related functions.
import sys  # It allows to operate on the interpreter as it provides access to variables and functions that interact strongly with the interpreter.

parser = argparse.ArgumentParser(
    prog="Variant_Calling_MUTEC2_ANNOVAR.py",
    description="""Variant calling - Mate pairs""",
)

parser.add_argument(
    "-o",
    "--outDir",
    dest="outDir",
    action="store",
    required=True,
    help="Path to out Working Directory.",
)

parser.add_argument(
    "-bam",
    "--bam",
    dest="bam",
    action="store",
    help="bam file generated in the Fastq_Bam_AC2.py",
)

options = parser.parse_args()

outDir = options.outDir
bam = options.bam

start_time = time.time()

sample = bam.split("/")[1].replace(".bam", "").split("_")[0]

LOG_FILE_OUT = open(outDir + "/log_VC/log_out_" + sample + ".txt", "w")

LOG_FILE_ERR = open(outDir + "/log_VC/log_err_" + sample + ".txt", "w")

bashArguments = (
    "java -Xmx16G -jar /apps/GATK/4.1.8.1/gatk-package-4.1.8.1-local.jar BaseRecalibrator -I "
    + bam
    + " -L Info_VC/BEDfile_Mutations_Panel_Design_no11q_168gens.bed -O Mutect2/"
    + sample
    + "_recal.txt -R Ref_Genome/genome.fa --known-sites Info_VC/dbsnp_138.hg19.vcf --known-sites Info_VC/1000G_phase1.indels.hg19.vcf"
)

e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

if e != 0:
    LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
    LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
    sys.exit(1)
bashArguments = (
    "java -Xmx16G -jar /apps/GATK/4.1.8.1/gatk-package-4.1.8.1-local.jar ApplyBQSR -I "
    + bam
    + " -L Info_VC/BEDfile_Mutations_Panel_Design.bed -O Mutect2/"
    + sample
    + "_BQSR.bam -R Ref_Genome/genome.fa --bqsr-recal-file Mutect2/"
    + sample
    + "_recal.txt"
)

if e != 0:
    LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
    LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
    sys.exit(1)
bashArguments = (
    "java -Xmx16G -jar /apps/GATK/4.1.8.1/gatk-package-4.1.8.1-local.jar Mutect2 -I Mutect2/"
    + sample
    + "_BQSR.bam -L Info_VC/BEDfile_Mutations_Panel_Design.bed -O Mutect2/"
    + sample
    + ".vcf -R Ref_Genome/genome.fa -tumor "
    + sample
    + " --native-pair-hmm-threads 8 --max-reads-per-alignment-start 0"
)

e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

if e != 0:
    LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
    LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
    sys.exit(1)
bashArguments = (
    "java -Xmx16G -jar /apps/GATK/4.1.8.1/gatk-package-4.1.8.1-local.jar FilterMutectCalls -V Mutect2/"
    + sample
    + ".vcf -O Mutect2/"
    + sample
    + "_filt.vcf --native-pair-hmm-threads 8 --reference Ref_Genome/genome.fa"
)

e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

if e != 0:
    LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
    LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
    sys.exit(1)

    LOG_FILE_OUT.write(
        "JOB DONE. Exectuion time in minutes: %s\n"
        % (round((time.time() - start_time) / 60, 4))
    )
    LOG_FILE_OUT.close()
    LOG_FILE_ERR.close()
