# DADA2-SA paper

Denis Mestivier - 2021

## Environnment creation

We used a conda environment:

```bash
conda create -n dada2
conda activate dada2
conda install -c bioconda bioconductor-dada2
conda install -c conda-forge r-ggplot2
conda install -c r r-tidyr 
conda install r-pheatmap
```

```bash
conda env export > dada2.yml
conda env create -f dada2.yml
```

## DADA

Raw data not provided here because :

```
REF:
   "Characterization of biliary microbiota dysbiosis in extrahepatic cholangiocarcinoma"
   Massa Saab, Denis Mestivier, Masoudreza Sohrabi, Christophe Rodriguez, Mahmood Reza Khonsari, Amirhossein Faraji, Iradj Sobhani
   PLoS One,  2021 Mar 9;16(3):e0247798.
   PMID: 33690612

The 16S rRNA amplicon sequencing data are available from the European
Nucleotide Archive (ENA) database (http://www.ebi.ac.uk/ena) under the
accession number PRJEB43183.
```

## Scripts

Different sets of parameters are addressed in this work.
Each set is tested in the `./PDS` directory.

Scripts for the figures are in the `./Figures` directory.

