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

	#L = which( colnames(mydf) == "Trunclen" )
	#colnames(mydf)[L] = "PRM_TRUNCLEN"

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

listOfPDS = list.files( path="../../PDS/Error-dada/", pattern="_PDS.csv", full.names=T )

# import data
m = c()

for( fn in listOfPDS )
{
	tmp = read.csv( fn, h=T, sep=",", stringsAsFactors=T )
	tmp = prep( tmp )
	m = rbind(m,tmp)
}

# wide to long
gathercols = c("pFiltered", "pMerged", "pNonchim" )
m.long = gather( m, stats, percent, gathercols ) 

# Selection subset

#L1 = which( m.long$PRM_TRUNCLEN %in% c("TruncLen_0-0", "TruncLen_250-250", "TruncLen_245-245", "TruncLen_245-235", "TruncLen_240-240", "TruncLen_235-235" ) )
#L2 = which( m.long$PRM_TRUNCLEN %in% c("TruncLen_0-0", "TruncLen_240-240", "TruncLen_245-235", "TruncLen_245-240", "TruncLen_240-230", "TruncLen_240-225" )  )

# Rename factors for visualisation
#m.long$Trunclen = m.long$PRM_TRUNCLEN
#levels(m.long$Trunclen) = gsub( "TruncLen_", "", levels(m.long$PRM_TRUNCLEN ) )

# Visualisation
ggplot( m.long, aes(x=nbases, y=percent, fill=stats ) ) + geom_boxplot() + ylim(0,100)

ggsave( "fig_truncLen-left.pdf", scale = 1, width = 21, height = 29.7, units = c("cm"  ),   dpi = 300, limitsize = TRUE )

