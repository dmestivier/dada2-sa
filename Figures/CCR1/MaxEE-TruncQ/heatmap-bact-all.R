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
pathFile   = "../../../PDS/CCR1/FilterTrim/MaxEE-TruncQ/"

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

	# sort levels
	#levels(pds$PRM_TRUNCQ) = c("1","2","4","8","16" )

	# Normalise into percent genus
	percent = apply( counts, 2, function(x){ 100.0*x / sum(x) } )

	# Filter non-genus taxa
	L       = which( is.na( taxa[,"Genus"] ) )
	if( length(L) > 0)
	{
		percent = percent[ -L, ]
		taxa    = taxa  [ -L, ]
	}


	rownames(percent) = taxa[ ,"Genus" ]
	# RColorBrewer

	mycol = colorRampPalette( brewer.pal(n = 7, name = "YlGnBu" ))

	# heatmap

	# cols annotation
	acols = pds[, c( "PRM_MAXEE", "PRM_TRUNCQ") ]
	#acols = pds[, c("PRM_TRUNCQ", "PRM_MAXEE") ]
	rownames(acols) = pds[,1]

	otfnALL  = paste( samplename, "_fig_maxEE-truncQ-heatmap-ALL.pdf", sep="" )
	otfnHIGH = paste( samplename, "_fig_maxEE-truncQ-heatmap-HIGH.pdf", sep="" )
	otfnLOW  = paste( samplename, "_fig_maxEE-truncQ-heatmap-LOW.pdf", sep="" )
	XX1=1
	XX2=min(20,nrow(percent))

	pheatmap( percent, cluster_cols=F, cluster_rows=F, gaps_col=c(5,10,15,20,25), angle_col=c("315"), annotation_col=acols, show_colnames=T, color=mycol(100), filename=otfnALL, width=8.25, height=11.75)
	
	pheatmap( percent[XX1:XX2,], cluster_cols=F, cluster_rows=F, gaps_col=c(5,10,15,20,25), angle_col=c("315"), annotation_col=acols, show_colnames=T, color=mycol(100), filename=otfnHIGH, width=8.25, height=11.75)

	if( nrow(percent)>10 )
	{	
		pheatmap( percent[(nrow(percent)-10):nrow(percent),], cluster_cols=F, cluster_rows=F, gaps_col=c(5,10,15,20,25), angle_col=c("315"), annotation_col=acols, show_colnames=T, color=mycol(100), filename=otfnLOW, width=8.25, height=11.75)
	}
}

