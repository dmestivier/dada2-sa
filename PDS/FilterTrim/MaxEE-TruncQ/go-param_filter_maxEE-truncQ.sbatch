#!/bin/bash
#SBATCH -p ALL 
#SBATCH -J slurm-dada2-error
#SBATCH -e slurm-dada2-error-%j
#SBATCH -o slurm-dada2-error-%j

#conda activate dada2

Rscript /work/dmestivier/Travail/ec2m3.paper.dada2/PDS/FilterTrim/MaxEE-TruncQ/param_filter_maxEE-truncQ.R $1
