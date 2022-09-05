# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 14:35:24 2022

@author: User
"""

import subprocess
import argparse
import time
import sys

parser = argparse.ArgumentParser(
    prog="Interval_list.py", description="""Interval list"""
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
    "-inp",
    "--inputfile",
    dest="inp",
    action="store",
    required=True,
    help="Path to the Bam for which we want to test the Quality.",
)

parser.add_argument(
    "-sd",
    "--sd",
    dest="sd",
    action="store",
    required=True,
    help="Path to the Bait Intervals in a interval_list format.",
)


options = parser.parse_args()

outDir = options.outDir
inp = options.inp
sd = options.sd

start_time = time.time()

LOG_FILE_OUT = open(outDir + "/log/log_out_il.txt", "w")

LOG_FILE_ERR = open(outDir + "/log/log_err_il.txt", "w")

bashArguments = (
    "java -Xmx16G -jar /apps/PICARD/2.24.0/bin/picard.jar BedToIntervalList  I="
    + inp
    + " O=interval_list_KIEL.interval_list SD="
    + sd
    + ""
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
