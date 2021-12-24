###
### Read counts/taxa & statistical output from Dada2 pipeline
### Produce heatmap 
###
### Mestivier Denis - November, 2021
###

##############################################
###
### Env
###
##############################################

rm(list=ls())
library(pheatmap)

##############################################
###
### main
###
##############################################

# filenames
samplename = "14SA_S77"
samplename = "19KA_S72"
fn.counts  = paste( samplename, "_counts.csv", sep="" )
fn.taxa    = paste( samplename, "_taxa.csv", sep="" )
fn.pds     = paste( samplename, "_PDS.csv", sep="" )

# import data
counts = read.csv( fn.counts, h=T, sep=",", stringsAsFactors=T, row.names=1 )
taxa   = read.csv( fn.taxa  , h=T, sep=",", stringsAsFactors=T, row.names=1 )
pds    = read.csv( fn.pds   , h=T, sep=",", stringsAsFactors=T )

# Normalise into percent genus
percent = apply( counts, 2, function(x){ 100.0*x / sum(x) } )

# Filter non-genus taxa
L       = which( is.na( taxa[,"Genus"] ) )
percent = percent[ -L, ]
taxa    = taxa  [ -L, ]

rownames(percent) = taxa[ ,"Genus" ]

# heatmap

pheatmap( percent, cluster_cols=F, cluster_rows=F, gaps_col=c(8,16,24,32), angle_col=c("315") )
