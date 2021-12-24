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

# list of sample
lof  = read.csv("list-of-files-ccr1", h=F, sep="\t" )
colnames(lof)[1] = c("sn")

# filename
pathFile   = "../../CCR1/FilterTrim/Trunclen/"

# pour chaque sample

for( samplename in lof[, 1] )
{
	fn.counts  = paste( pathFile, samplename, "_counts.csv", sep="" )
	fn.taxa    = paste( pathFile, samplename, "_taxa.csv", sep="" )
	fn.pds     = paste( pathFile, samplename, "_PDS.csv", sep="" )

	# import data
	counts = read.csv( fn.counts, h=T, sep=",", stringsAsFactors=T, row.names=1 )
	taxa   = read.csv( fn.taxa  , h=T, sep=",", stringsAsFactors=T, row.names=1 )
	pds    = read.csv( fn.pds   , h=T, sep=",", stringsAsFactors=T )

	# Rename factors for visualisation
	levels(pds$Trunclen) = gsub( "TruncLen_", "", levels(pds$Trunclen ) )

	# Selection subset
	L1      = which( pds$Trunclen %in% c("0-0", "250-250", "245-245", "245-235", "240-240", "235-235" ) )
	pds     = pds[ L1, ]
	counts  = counts[ , L1 ]

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

	THR = 0.001
	p2 = percent
	p2[ which( p2 < THR ) ] = NA

	percent = p2

	# cols annotation
	acols = as.data.frame( pds[, c("Trunclen") ] )
	rownames(acols) = pds[,1]
	colnames(acols) = c("truncLen")

	# 
	XX1=1
	XX2=15
	XX3=nrow(percent) - XX2

	otfnALL  = paste( samplename, "_fig_truncLen-heatmap-full.pdf", sep="" )
	otfnHIGH = paste( samplename, "_fig_truncLen-heatmap-high.pdf", sep="" )
	otfnLOW  = paste( samplename, "_fig_truncLen-heatmap-low.pdf", sep="" )

	pheatmap( percent, cluster_cols=F, cluster_rows=F, gaps_col=c(5), angle_col=c("315"), annotation_col=acols, show_colnames=F, color=mycol(100), filename=otfnALL, width=8.25, height=11.75, na_col="red" )

	pheatmap( percent[XX1:XX2,], cluster_cols=F, cluster_rows=F, gaps_col=c(5), angle_col=c("315"), annotation_col=acols, show_colnames=F, color=mycol(100), filename=otfnHIGH, width=8.25, height=11.75, na_col="red" )

	pheatmap( percent[XX3:nrow(percent),], cluster_cols=F, cluster_rows=F, gaps_col=c(5), angle_col=c("315"), annotation_col=acols, show_colnames=F, color=mycol(100), filename=otfnLOW, width=8.25, height=11.75, na_col="red" )
}
