#!/bin/bash
#SBATCH --job-name="VARIATN CALLER MUTECT2"
#SBATCH -D .
#SBATCH --output=Mutect2.out
#SBATCH --error=Mutect2.err
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --qos=debug

module load greasy/2.2
module load java/10.0.1
module load gatk/4.1.8.1
module load perl/5.26.2

greasy 2_MUTECT2_FINAL.txt
