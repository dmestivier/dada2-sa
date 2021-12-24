library(dada2)
library(dplyr)

rm( list=ls() )

path<-"."

#pathDB <- "/home/dmestivier/Publications/2020-2021/ec2m3.paper.dada2/DB"
#pathFQ <- "/home/dmestivier/Publications/2020-2021/ec2m3.paper.dada2/Data"
#pathRES<- "/home/dmestivier/Publications/2020-2021/ec2m3.paper.dada2/PDS/Error-dada"
pathDB <- "/work/dmestivier/Travail/ec2m3.paper.dada2/DB"
pathFQ <- "/work/dmestivier/Travail/ec2m3.paper.dada2/Data"
pathRES<- "/work/dmestivier/Travail/ec2m3.paper.dada2/PDS/Error-dada"

##NOM FICHIER 
T1<-Sys.time()

sample.name <- "19KA_S72"
fnFs <- paste( pathFQ, "/", sample.name, "_L001_R1_001.fastq", sep="" )
fnRs <- paste( pathFQ, "/", sample.name, "_L001_R2_001.fastq", sep="" )

param<-NULL

nbases<-c( 1e4, 1e6, 1e8, 1e10)
truncL<-c( c(245,235) )
MAX_CONSIST=c(5,10,20)               #learnError()
OMEGA_A=c(1e-40,1e-30,1e-20, 1e-10)  #dada()

for(i in 1:(length(truncL)/2))
{
	for(j in 1:length(nbases))
	{
		for( l in 1:length(MAX_CONSIST) )
		{
			for( k in 1:length(OMEGA_A) )
			{
				param<-rbind(param,list(c(truncL[2*(i-1)+1],truncL[2*(i-1)+2]),nbases[j], MAX_CONSIST[l], OMEGA_A[k]))
			}
		}
    }
}

param<-as.data.frame(param)
colnames(param)<-c("Tcl","nbases", "MaxConsist", "OmegaA" )

#######################################
table_taxa<-NULL

for (i in 1:nrow(param))
{
	cat("Run = ", i, "\n" )
	cat("\n")

	####pipeline dada ########################################

	filtFs <- file.path(pathFQ, "filtered2", paste0(sample.name, "_F_filt.fastq.gz"))
	filtRs <- file.path(pathFQ, "filtered2", paste0(sample.name, "_R_filt.fastq.gz"))

	out<-filterAndTrim( fnFs, filtFs, fnRs, filtRs, truncLen=param$Tcl[[i]],maxN=0, maxEE=c(2,2), 
                        truncQ=2, rm.phix=T,compress=TRUE, multithread = T )

	errF <- learnErrors( filtFs, multithread=TRUE, nbases=param$nbases[[i]], MAX_CONSIST=param$MaxConsist[[i]]  )
  	errR <- learnErrors( filtRs, multithread=TRUE, nbases=param$nbases[[i]], MAX_CONSIST=param$MaxConsist[[i]]  )
  
	#plotErrors(errF, nominalQ=TRUE)
  
	dadaFs <- dada( filtFs, err=errF, multithread=TRUE, OMEGA_A=param$OmegaA[[i]] )
  	dadaRs <- dada( filtRs, err=errR, multithread=TRUE, OMEGA_A=param$OmegaA[[i]] )

  	mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, trimOverhang=F, verbose=TRUE )
  	seqtab <- makeSequenceTable(mergers)
  	print(dim(seqtab))
  	seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  	taxa<- assignTaxonomy(seqtab.nochim, paste( pathDB, "/silva_nr_v132_train_set.fa.gz", sep="" ), multithread=T,minBoot =50,outputBootstraps = F, tryRC=F )
  	#taxa <- addSpecies(taxa, "../../../silva_species_assignment_v132.fa.gz")
  
	##############################
 
  	##Lecture de la table d'assignation taxonomique 
  	taxa.print <- taxa 
  	rownames(taxa.print) <- NULL
  	taxa_data<-as.data.frame(taxa.print)
  	nochim<-seqtab.nochim
  	colnames(seqtab.nochim)<-NULL
  	ind<-which(is.na(taxa.print[,2]),arr.ind=F)
  	ind_e<-which(taxa.print[,1]=="Eukaryota")
  	taxa_na<-taxa.print
  	taxa_na<-taxa_na[-c(ind,ind_e),]
  	##concaténation des colonnes , chemins taxonomiques 
  	chemin_tax<-c()
  	nom_col<-c()
  	for(j in 1:nrow(taxa_data))
	{
    	for(k in 2:7)
		{
      		nom_col<-paste(nom_col,taxa_data[j,k],sep="")
    	}
    	chemin_tax <-c(chemin_tax,nom_col)
    	nom_col<-c()
  	}
  	##ajout des colonnes chemin tax et du nombre de reads par asv 
  	taxa_data<-cbind(taxa_data,chemin_tax,seqtab.nochim[1,])
  	colnames(taxa_data)[ncol(taxa_data)]<-paste("track",i,sep="")
  	colnames(taxa_data)[(ncol(taxa_data)-1)]<-c("nom")
  	if(i>1)
	{
    	##taxa
    
    	#fusion des tables track i-1 et track i 
    	table_fusion<-bind_rows(table_taxa,taxa_data)
    	##table comptabilisé les reads pour un chemin taxonomique
    	table_count<-NULL  
    	##table taxonomique à construire
    	table_tx<-NULL  
    	name<-unique(table_fusion$nom)
    	vec_sum<-c()
    	for(j in 1:length(name))
		{
      		## recuperation de toutes les lignes ayant le meme chemin taxonomique
      		table_count<-rbind(table_count,filter(table_fusion,nom==name[j]))
      		## on recupere les colonnes de la table taxonomique
      		colonne_ta<-table_count[1,1:7]
      		## on somme les colonnes de reads de chaque track
      		for(k in 8:ncol(table_fusion))
			{
        		colonne_ta<-cbind(colonne_ta,sum(table_count[,k],na.rm=T))
      		}
      		##on ajoute la nouvelle ligne du chemin taxonomique avec les reads sommés
      		table_tx<-rbind(table_tx,colonne_ta)
      		## on réinitialise l'objet table_count 
      		table_count<-NULL
    	}
    
		table_taxa<-table_tx
    	##nom de colonnes
    	for(k in 8:ncol(table_fusion))
		{
      		colnames(table_taxa)[k]<-paste("track",(k-7),sep="")
    	}
  	}else{ 
    	##chemin taxonomique
    	table_taxa<-taxa_data
  	}
  	rownames(out)<-NULL
  	na<-nrow(taxa_na)
  
	if(is.null(nrow(taxa_na)))
	{
    	na<-0
  	}
  	
	track <- cbind(paste("track",i,sep=""),out,sum(getUniques(dadaFs)),sum(getUniques(dadaRs)),sum(getUniques(mergers)), rowSums(seqtab.nochim),na,length(unique(c(ind_e,ind))),
                 paste("TruncLen_",paste(param$Tcl[[i]][1],param$Tcl[[i]][2],sep="-"),sep=""),paste("nbases_",param$nbases[[i]],sep=""),paste("MaxConsist_",param$MaxConsist[[i]], sep=""), paste( "OmegaA_", param$OmegaA[[i]] , sep="") )
  
	colnames(track) <- c("","input","filtered", "denoisedF", "denoisedR", "merged", "nonchim","assignes","badassignes","Trunclen","nbases","MaxConsist", "OmegaA")
  
	##changement du nom des colonnes de 
	for (j in 1:ncol(nochim))
	{
    	colnames(nochim)[j]<-c(paste("asv",j,sep=""))
  	}
  
	##ecriture dans des fichiers 
  	if(i>1)
	{
		write.table(track,paste(sample.name,"_error-dada_PDS.csv",sep=""),append=T,sep=",",qmethod="double",row.names = F,col.names = F)
  	}else{
    	write.table(track,paste(sample.name,"_error-dada_PDS.csv",sep=""),append=F,sep=",",qmethod="double",row.names = F,col.names = T)
  	}
}

asv<-c()  
for(j in 1:nrow(table_taxa))
{
	asv<-c(asv,paste("asv",j,sep=""))
}

table_taxa<-cbind(asv,table_taxa)
colnames(table_taxa)[1]<-c("asv")

write.table(table_taxa[,1:7],paste(sample.name,"_error-dada_taxa.csv",sep=""),append=F,sep=",",qmethod="double",row.names = F,col.names = T)

t<-cbind(table_taxa[,1],table_taxa[,9:ncol(table_taxa)])

write.table(t,paste(sample.name,"_error-dada_counts.csv",sep=""),append=F,sep=",",qmethod="double",row.names = F,col.names = T)

T2<-Sys.time()
Tdiff<-difftime(T2, T1)

cat("End of pipeline\n" )
