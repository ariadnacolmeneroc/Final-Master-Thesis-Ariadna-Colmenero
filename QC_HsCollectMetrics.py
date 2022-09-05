# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 14:34:55 2022

@author: User
"""

import subprocess
import argparse
import time
import sys

parser = argparse.ArgumentParser(
    prog="QC_HsCollect.py", description="""BAM FILES QUALITY CONTROL"""
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
    "--bamfile",
    dest="bam",
    action="store",
    required=True,
    help="Path to the Bam for which we want to test the Quality.",
)

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
)


options = parser.parse_args()

outDir = options.outDir
bam = options.bam
ref = options.ref
bait = options.bait

start_time = time.time()

sample = bam.split("/")[-1].split("_")[0]

LOG_FILE_OUT = open(outDir + "/log/log_out_" + sample + ".txt", "w")

LOG_FILE_ERR = open(outDir + "/log/log_err_" + sample + ".txt", "w")

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
)

e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

if e != 0:
    LOG_FILE_ERR.write("ERROR IN COMMAND: " + bashArguments + "\n")
    LOG_FILE_OUT.write("ERROR IN COMMAND: " + bashArguments + "\n")
    sys.exit(1)
