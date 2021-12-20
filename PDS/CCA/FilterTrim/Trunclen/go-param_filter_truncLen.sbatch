#!/bin/bash
#SBATCH -p ALL 
#SBATCH -J slurm-dada2-error
#SBATCH -e slurm-dada2-error-%j
#SBATCH -o slurm-dada2-error-%j

#conda activate dada2

Rscript /work/dmestivier/Travail/ec2m3.paper.dada2/PDS/FilterTrim/Trunclen/param_filter_truncLen.R $1
