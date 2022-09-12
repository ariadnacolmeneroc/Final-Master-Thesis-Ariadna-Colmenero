
############################### SH SCRIPT INTERVAL LIST ########################################

# Ariadna Colmenero Cobo de Guzman
# Master's Final Thesis
# 2021-2022

# A.sh file is a scripting language commands file that contains computer programmes that can be run by the Unix shell. This is one of the three files we need to be able to run Interval_list.py:

#!/bin/bash # This command is telling our system to use bash as the default shell.
#SBATCH --job-name="FASTQTOBAM" # Job name for the commit. It has to be descriptive so we can know which is each job in the queue.
#SBATCH --output=FASTQBAM.out # The output will be saved as it is produced. We can open this file and follow the job.
#SBATCH --error=FASTQBAM.err # If there is an error, it will tell us exactly what the error is and on which line it is located.
#SBATCH --ntasks=15 # Number of tasks performed in parallel. For example, in this case the .PY is applied for the formation of two different BAMs.
#SBATCH --time=2:00:00 # Time we give for the script to run. If it reaches the time limit, it will stop. 
                        # Further consideration: If it stops and the job has not yet finished, look at the OUTPUT file. We must check that it has been done completely. 
#SBATCH --cpus-per-task=8 # It will ensure that the task gets allocated to the same node.
#SBATCH --qos=debug #  Quality of Service (QOS). As we are launching a 2 h job, we can batch it in the fast queue.

# Next, we need to load all the modules that are needed to run the corresponding PY file:

module load greasy/2.2
module load java/10.0.1
module load picard/2.24.0
module load R/3.5.0
module load pandoc

greasy interval_list.txt # We use it to spawn tasks to remote nodes inside a job. So, this job is going to use this TXT request to locate the required files via the paths it contains.
