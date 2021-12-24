###
### Read statiscal output from Dada2 pipeline
### Meerge samples
### Produce boxplot
###
###   Rscript draw_percentage.R
###
### Mestivier Denis - November, 2021
###

##############################################
###
### Env
###
##############################################

rm(list=ls())
library(ggplot2)
library(tidyr)
library(reshape2)

#############################################
###
### Prepare dataframe
###
#############################################

prep = function( mydf )
{
	colnames(mydf)[1] = "track"

	# Some factor should be reorderd / renamed
	mydf$PRM_TRUNCQ = factor( mydf$PRM_TRUNCQ, levels=c("Truncq_1", "Truncq_2", "Truncq_4", "Truncq_8", "Truncq_16" ) )

	L = which( colnames(mydf) == "Trunclen" )
	colnames(mydf)[L] = "PRM_TRUNCLEN"

	# Normalisation
	mydf$pFiltered  = 100.0 * mydf$filtered  / mydf$input
	mydf$pDenoisedF = 100.0 * mydf$denoisedF / mydf$input
	mydf$pDenoisedR = 100.0 * mydf$denoisedR / mydf$input
	mydf$pMerged    = 100.0 * mydf$merged    / mydf$input
	mydf$pNonchim   = 100.0 * mydf$nonchim   / mydf$input
	
	mydf
}

##############################################
###
### main
###
##############################################

# list of filenames
listOfPDS = list.files( path="../../ITMO16S/FilterTrim/MaxEE-TruncQ", pattern="_PDS.csv", full.names=T )

# import data
m = c()

for( fn in listOfPDS )
{
	tmp = read.csv( fn, h=T, sep=",", stringsAsFactors=T )
	tmp = prep( tmp )
	m   = rbind(m,tmp)
}

# Visualisation
ggplot( m, aes( x=PRM_TRUNCQ, y=pNonchim, fill=PRM_MAXEE)) + geom_boxplot() + ylim(0,100)

ggsave( "fig_maxEE-truncQ-left.pdf", scale = 1, width = 21, height = 29.7, units = c("cm"  ),   dpi = 300, limitsize = TRUE )

