###
### Generate figure for read quality profiles
###
### Mestivier Dens - Novembre, 2021
###

# clean env
rm( list=ls() )

# libraries needed
library(dada2)
library(ggplot2)

# path

#pathFQ <- "/work/dmestivier/Travail/ec2m3.paper.dada2/Data"
pathFQ <- "/home/dmestivier/Publication/2021-2022/ec2m3.paper.dada2/Data"

## for the whole cohort 
#fnFs <- sort(list.files(pathFQ, pattern="_L001_R1_001.fastq", full.names = TRUE))
#fnRs <- sort(list.files(pathFQ, pattern="_L001_R2_001.fastq", full.names = TRUE))

## Only one file for publication

fnFs <- "/home/dmestivier/Publication/2021-2022/ec2m3.paper.dada2/Data/19KA_S72_L001_R1_001.fastq"
fnRs <- "/home/dmestivier/Publication/2021-2022/ec2m3.paper.dada2/Data/19KA_S72_L001_R2_001.fastq"

### Inspect read quality profiles

plotQualityProfile(fnFs[1:2])
ggsave( "fig_readQualProfDADA2_R1.pdf", plot = last_plot(), scale = 1, width = NA, height = NA, units = c("cm"), dpi = 300, limitsize = TRUE )

plotQualityProfile(fnRs[1:2])
ggsave( "fig_readQualProfDADA2_R2.pdf", plot = last_plot(), scale = 1, width = NA, height = NA, units = c("cm"), dpi = 300, limitsize = TRUE )



