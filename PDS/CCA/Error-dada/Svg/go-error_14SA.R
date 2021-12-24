#!/bin/bash
#SBATCH -p ALL 
#SBATCH -J dada2-error
#SBATCH -e dada2-error-%j
#SBATCH -o dada2-error-%j

#conda activate dada2

Rscript /work/dmestivier/Travail/ec2m3.paper.dada2/PDS/Error/error_14SA.R
