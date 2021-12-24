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
	L = which( colnames(mydf) == "Trunclen" )
	colnames(mydf)[L] = "PRM_TRUNCLEN"
	
	# reorder factor
	mydf$PRM_TRUNCLEN = factor( mydf$PRM_TRUNCLEN, 
		levels=
			   c( 
				 "TruncLen_0-0",
				"TruncLen_250-250",
				"TruncLen_245-245",
				"TruncLen_240-240",
				"TruncLen_235-235",

				"TruncLen_235-220" ,
				"TruncLen_235-225",
				"TruncLen_235-230",

				"TruncLen_240-220",
				"TruncLen_240-225",
				"TruncLen_240-230",
				"TruncLen_240-235",

				"TruncLen_245-225",
				"TruncLen_245-230",
				"TruncLen_245-235",
				"TruncLen_245-240",

				"TruncLen_250-225",
				"TruncLen_250-230",
				"TruncLen_250-235",
				"TruncLen_250-240",
				"TruncLen_250-245" ) )


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
f3 = "24MA_S64_PDS.csv"

# import data
m1 = read.csv( f1, h=T, sep=",", stringsAsFactors=T )
m2 = read.csv( f2, h=T, sep=",", stringsAsFactors=T )
m3 = read.csv( f2, h=T, sep=",", stringsAsFactors=T )

m1 = prep( m1 )
m2 = prep( m2 )
m3 = prep( m3 )

# Agregation
mm = rbind(m1,m2)
mm = rbind(mm,m3)

# wide to long
gathercols = c("pFiltered", "pMerged", "pNonchim" )
#gathercols = c("pFiltered", "pDenoisedF", "pDenoisedR", "pMerged", "pNonchim" )
mm.long = gather( mm, stats, percent, gathercols ) 

# Filter
L = which( mm$PRM_TRUNCLEN %in% c(
							"TruncLen_0-0",
							"TruncLen_250-250",
							"TruncLen_245-245",
							"TruncLen_240-240",
							"TruncLen_235-235" ) )

mm.ok = mm[ L, ]
mm.ok.long = gather( mm.ok, stats, percent, gathercols )

mm.off = mm[ -L, ]
mm.off.long = gather( mm.off, stats, percent, gathercols )

# Visualisation
ggplot( mm.ok, aes(x=PRM_TRUNCLEN, y=pFiltered ) ) + geom_boxplot()
