#!/bin/bash
#SBATCH -p ALL 
#SBATCH -J slurm-dada2-filt-trim
#SBATCH -e slurm-dada2-filt-trim-%j
#SBATCH -o slurm-dada2-filt-trim-%j

#conda activate dada2

Rscript /work/dmestivier/Travail/ec2m3.paper.dada2/PDS/FilterTrim/param_filter_14SA.R
