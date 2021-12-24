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


REF:
   "The Limits and Avoidance of Biases in Metagenomic Analyses of Human Fecal Microbiota"
   Emma Bergsten, Denis Mestivier, Iradj Sobhani
   Microorganisms, 2020 Dec 9;8(12):1954.
   doi: 10.3390/microorganisms8121954.
   PMID: 33317070

The CCR1/16S were downloaded from the European Nucleotide Archive (ENA) database
(http://www.ebi.ac.uk/ena) under the access number ERP005534.
```

ITMO/16S data are not yet published.

Data are gzipped but should be un-gzipped to reproduce the figures.

## Scripts

Different sets of parameters are addressed in this work.
Each parameter set tested is in the `./PDS` directory:

- MaxEE and truncQ
- truncLen
- LearnError and dada
- assignTaxonomy

Scripts for the figures are in the `./Figures` directory.

- ReadQualityProfiles
- MaxEE-TruncQ
- TruncLen
- Error-Dada
- assigTaxo

