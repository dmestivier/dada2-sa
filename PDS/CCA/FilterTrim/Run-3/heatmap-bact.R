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
library(RColorBrewer)

##############################################
###
### main
###
##############################################

# filenames
samplename = "19KA_S72"
samplename = "14SA_S77"
#samplename = "24MA_S64"

fn.counts  = paste( samplename, "_counts.csv", sep="" )
fn.taxa    = paste( samplename, "_taxa.csv", sep="" )
fn.pds     = paste( samplename, "_PDS.csv", sep="" )

# import data
counts = read.csv( fn.counts, h=T, sep=",", stringsAsFactors=T, row.names=1 )
taxa   = read.csv( fn.taxa  , h=T, sep=",", stringsAsFactors=T, row.names=1 )
pds    = read.csv( fn.pds   , h=T, sep=",", stringsAsFactors=T )

# sort levels
#levels(pds$PRM_TRUNCQ) = c("1","2","4","8","16" )

# Normalise into percent genus
percent = apply( counts, 2, function(x){ 100.0*x / sum(x) } )

# Filter non-genus taxa
L       = which( is.na( taxa[,"Genus"] ) )
percent = percent[ -L, ]
counts  = counts [ -L, ]
taxa    = taxa  [ -L, ]

rownames(percent) = taxa[ ,"Genus" ]

# Filter truncLen values

L = which( pds$Trunclen %in% c(
	"TruncLen_0-0",
	"TruncLen_250-250",
	"TruncLen_245-245",
	"TruncLen_240-240",
	"TruncLen_235-235" ) )

percent = percent[ , L ]
counts  = counts [ , L ]
pds     = pds    [ L, ]

# Genus at a very very low lever (<0.01%) are 
# converted to NA (for pheatmap/na_col="red")

Z = which( percent < 0.01, arr.ind=T )
percent[Z] = NA

# RColorBrewer
mycol = colorRampPalette( brewer.pal(n = 7, name = "YlGnBu" ))

# heatmap

# cols annotation
acols = as.data.frame( pds[, c("Trunclen") ] )
rownames(acols) = pds[,1]
colnames(acols) = c("truncLen")

# ou prend le log2(c + 1 )
pheatmap( percent, cluster_cols=F, cluster_rows=F, angle_col=c("315"), annotation_col=acols, show_colnames=F, color=mycol(100), na_col="red" )

XX=10
#pheatmap( percent[1:XX,], cluster_cols=F, cluster_rows=F, gaps_col=c(5,10,15,20,25), angle_col=c("315"), annotation_col=acols, show_colnames=F, color=mycol(100))
#pheatmap( percent[10:nrow(percent),], cluster_cols=F, cluster_rows=F, gaps_col=c(5,10,15,20,25), angle_col=c("315"), annotation_col=acols, show_colnames=F, color = mycol(100) ) 
