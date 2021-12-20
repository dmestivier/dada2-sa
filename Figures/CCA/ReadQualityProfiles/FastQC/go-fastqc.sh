#!/bin/bash

#conda activate dada2

for fn in ../../../Data/*.fastq
do
	fastqc $fn -o .
done
