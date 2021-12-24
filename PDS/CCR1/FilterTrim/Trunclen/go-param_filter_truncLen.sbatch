#!/bin/bash
#SBATCH -p ALL 
#SBATCH -J slurm-ccr1-trunclen
#SBATCH -e slurm-ccr1-trunclen-%j
#SBATCH -o slurm-ccr1-trunclen-%j

#conda activate dada2

Rscript /work/dmestivier/Travail/ec2m3.paper.dada2/CCR1/FilterTrim/Trunclen/param_filter_truncLen-FASTA.R $1
