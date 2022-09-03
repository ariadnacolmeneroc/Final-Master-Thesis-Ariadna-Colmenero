#!/bin/bash
#SBATCH --job-name="FASTQTOBAM"
#SBATCH -D .
#SBATCH --output=FASTQBAM.out
#SBATCH --error=FASTQBAM.err
#SBATCH --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
######SBATCH --qos=debug

module load greasy/2.2
module load java/10.0.1
module load samtools/1.9
module load bwa/0.7.17
module load picard/2.24.0

greasy FASTQTOBAM.txt


