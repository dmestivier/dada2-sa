###
### Read statiscal output from Dada2 pipeline
### Focus on bacteria genus
###
### Produce boxplot
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
	mydf$PRM_TRUNCQ = factor( mydf$PRM_TRUNCQ, levels=c("Truncq_2", "Truncq_4", "Truncq_7", "Truncq_11", "Truncq_15" ) )

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

# filename
fn.pds = "19KA_S72_PDS.csv"
fn.counts = "19KA_S72_counts.csv"
fn.taxa = "19KA_S72_taxa.csv" 

# import data
tg = read.csv( fn.pds, h=T, sep=",", stringsAsFactors=T )
counts = read.csv( fn.counts, h=T, sep=",", stringsAsFactors=T, row.names=1 )
taxa   = read.csv( fn.taxa  , h=T, sep=",", stringsAsFactors=F, row.names=1 )

#
tg = prep( tg)

# wide to long for targetfile
#gathercols = c("pFiltered", "pDenoisedF", "pDenoisedR", "pMerged", "pNonchim" )
#mm.long = gather( mm, stats, percent, gathercols ) 

# Merge counts+taxo to targetfile
tcounts = t(counts)
colnames(tcounts)=taxa[, "Genus"] 
tg = cbind( tg, tcounts )

# remove "NA" Genus
L = which( is.na( colnames(tg)) )
tg = tg[ , -L ]

# Visualisation
#ggplot( m1, aes(x=PRM_MAXEE, y=pFiltered, fill=PRM_TRUNCQ ) ) + geom_boxplot()
ggplot( m1, aes(x=PRM_MAXEE, y=pFiltered ) ) + geom_boxplot()
#ggplot( m1, aes(x=PRM_MAXEE, y=pFiltered, fill=PRM_TRUNCQ ) ) + geom_boxplot() + facet_wrap( ~ PRM_TRUNCQ )
#ggplot( m1, aes(x=PRM_MAXEE, y=pFiltered, fill=PRM_TRUNCQ ) ) + geom_boxplot() + scale_y_continuous( name="Percent" )  + scale_x_discrete( name="MAX_EE", labels=c("A","B","C","D", "E" ) )

