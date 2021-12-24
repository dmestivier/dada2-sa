###
### Read statiscal output from Dada2 pipeline
### Meerge samples
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


# filenames
f1 = "19KA_S72_PDS.csv"
f2 = "14SA_S77_PDS.csv"


# import data
m1 = read.csv( f1, h=T, sep=",", stringsAsFactors=T )
m2 = read.csv( f2, h=T, sep=",", stringsAsFactors=T )

m1 = prep( m1 )
m2 = prep( m2 )

# Agregation
mm = rbind(m1,m2)

# wide to long
gathercols = c("pFiltered", "pDenoisedF", "pDenoisedR", "pMerged", "pNonchim" )
mm.long = gather( mm, stats, percent, gathercols ) 

# Visualisation
#ggplot( m1, aes(x=PRM_MAXEE, y=pFiltered, fill=PRM_TRUNCQ ) ) + geom_boxplot()
ggplot( m1, aes(x=PRM_MAXEE, y=pFiltered ) ) + geom_boxplot()
#ggplot( m1, aes(x=PRM_MAXEE, y=pFiltered, fill=PRM_TRUNCQ ) ) + geom_boxplot() + facet_wrap( ~ PRM_TRUNCQ )
#ggplot( m1, aes(x=PRM_MAXEE, y=pFiltered, fill=PRM_TRUNCQ ) ) + geom_boxplot() + scale_y_continuous( name="Percent" )  + scale_x_discrete( name="MAX_EE", labels=c("A","B","C","D", "E" ) )

