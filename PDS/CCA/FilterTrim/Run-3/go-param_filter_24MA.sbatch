#!/bin/bash
#SBATCH -p ALL 
#SBATCH -J dada2-filt-trim
#SBATCH -e dada2-filt-trim-%j
#SBATCH -o dada2-filt-trim-%j

conda activate dada2

Rscript /work/dmestivier/Travail/ec2m3.paper.dada2/PDS/FilterTrim/Run-3/param_filter_24MA.R
