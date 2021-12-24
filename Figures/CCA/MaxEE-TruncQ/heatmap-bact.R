###
### Read counts/taxa & statistical output from Dada2 pipeline
### Produce heatmap 
###
###  Rscript heatmap-bact.R
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
library(RColorBrewer)

##############################################
###
### main
###
##############################################

# filename

pathFile   = "../../../PDS/CCA/FilterTrim/MaxEE-TruncQ/"
samplename = "19KA_S72"
#samplename = "14SA_S77"
#samplename = "24MA_S64"

fn.counts  = paste( pathFile, samplename, "_counts.csv", sep="" )
fn.taxa    = paste( pathFile, samplename, "_taxa.csv", sep="" )
fn.pds     = paste( pathFile, samplename, "_PDS.csv", sep="" )

# import data
counts = read.csv( fn.counts, h=T, sep=",", stringsAsFactors=T, row.names=1 )
taxa   = read.csv( fn.taxa  , h=T, sep=",", stringsAsFactors=T, row.names=1 )
pds    = read.csv( fn.pds   , h=T, sep=",", stringsAsFactors=T )

# Rename factors for visualisation
#levels(pds$Trunclen) = gsub( "TruncLen_", "", levels(pds$Trunclen ) )

# Selection subset
#L1 = which( pds$Trunclen %in% c("0-0", "250-250", "245-245", "245-235", "240-240", "235-235" ) )
#L1 = which( pds$Trunclen %in% c("TruncLen_0-0", "TruncLen_250-250", "TruncLen_245-245", "TruncLen_245-235", "TruncLen_240-240", "TruncLen_235-235" ) )

#pds = pds[ L1, ]
#counts = counts[ , L1 ]

# sort levels
#levels(pds$PRM_TRUNCQ) = c("1","2","4","8","16" )

# Normalise into percent genus
percent = apply( counts, 2, function(x){ 100.0*x / sum(x) } )

# Filter non-genus taxa
L       = which( is.na( taxa[,"Genus"] ) )
percent = percent[ -L, ]
taxa    = taxa  [ -L, ]

rownames(percent) = taxa[ ,"Genus" ]

# RColorBrewer
mycol = colorRampPalette( brewer.pal(n = 7, name = "YlGnBu" ))

# heatmap

# genus abundance below THR is map as NA

THR = 0.01 
p2 = percent
p2[ which( p2 < THR ) ] = NA

# cols annotation
acols = as.data.frame( pds[, c("PRM_MAXEE", "PRM_TRUNCQ") ] )
rownames(acols) = pds[,1]

XX1=10

pheatmap( percent, cluster_cols=F, cluster_rows=F, angle_col=c("315"), annotation_col=acols, show_colnames=F, color=mycol(100), filename="fig_MaxEE-TruncQ-heatmap-full.pdf", width=8.25, height=11.75, na_col="black" )

pheatmap( percent[1:XX1,], cluster_cols=F, cluster_rows=F, angle_col=c("315"), annotation_col=acols, show_colnames=F, color=mycol(100), filename="fig_MaxEE-TruncQ-heatmap-high.pdf", width=8.25, height=11.75, na_col="black" )

#pheatmap( percent[10:30,], cluster_cols=F, cluster_rows=F, angle_col=c("315"), annotation_col=acols, show_colnames=F, color=mycol(100), filename="fig_MaxEE-TruncQ-dada-heatmap-medium.pdf", width=8.25, height=11.75, na_col="black" )

pheatmap( percent[(nrow(percent)-XX1):nrow(percent),], cluster_cols=F, cluster_rows=F, angle_col=c("315"), annotation_col=acols, show_colnames=F, color=mycol(100), filename="fig_MaxEE-TruncQ-heatmap-low.pdf", width=8.25, height=11.75, na_col="black" )


XX=10
#pheatmap( percent[1:XX,], cluster_cols=F, cluster_rows=F, gaps_col=c(5,10,15,20,25), angle_col=c("315"), annotation_col=acols, show_colnames=F, color=mycol(100))
#pheatmap( percent[10:nrow(percent),], cluster_cols=F, cluster_rows=F, gaps_col=c(5,10,15,20,25), angle_col=c("315"), annotation_col=acols, show_colnames=F, color = mycol(100) ) 
