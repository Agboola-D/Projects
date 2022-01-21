#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=156:59:59
#SBATCH --partition=smem
#SBATCH --ntasks=48
#SBATCH --job-name=Simulation2
#SBATCH --output=Data-Simulation2.%j.txt

module purge
source /curc/sw/anaconda3/latest
conda activate davidenv

Rscript Simulation2.R
