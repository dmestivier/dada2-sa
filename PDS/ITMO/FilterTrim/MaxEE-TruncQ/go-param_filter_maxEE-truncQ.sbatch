#!/bin/bash
#SBATCH -p ALL 
#SBATCH -J slurm-itmco-truncQ
#SBATCH -e slurm-itmco-truncQ-%j
#SBATCH -o slurm-itmoc-truncQ-%j

#conda activate dada2

Rscript /work/dmestivier/Travail/ec2m3.paper.dada2/ITMO16S/FilterTrim/MaxEE-TruncQ/param_filter_maxEE-truncQ.R $1

