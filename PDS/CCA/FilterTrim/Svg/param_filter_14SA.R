library(dada2)
library(dplyr)

rm( list=ls() )

path<-"."
pathFQ <- "/work/dmestivier/Travail/ec2m3.paper.dada2/Data"
pathDB <- "/work/dmestivier/Travail/ec2m3.paper.dada2/DB"
pathRES <-"/work/dmestivier/Travail/ec2m3.paper.dada2/PDS/FilterTrim"

##NOM FICHIER 
T1<-Sys.time()

fnFs <- paste( pathFQ, "/14SA_S77_L001_R1_001.fastq", sep="" )
fnRs <- paste( pathFQ, "/14SA_S77_L001_R2_001.fastq", sep="" )
sample.names <- "14SA_S77"

##PARAMETRES : FilterAndTRim 
##  maxEE,truncQ,rmPhix, TruncLen 

param<-NULL
maxEE<-c( c(1,1), c(2,2), c(3,3), c(2,5), c(3,5))
rm_phix<-c(T,F)
truncL<-c(c(0,0),c(250,230),c(230,210), c(210,190), c(180,160))
truncQ<-c(2,4,7,11,15)

for( i in 1:(length(maxEE)/2) )
{
  for( j in 1:length(rm_phix) )
  {
    for( k in 1:(length(truncL)/2) )
	{
      for( l in 1:length(truncQ) )
	  {
        param<-rbind(param,list(c(maxEE[2*(i-1)+1],maxEE[2*(i-1)+2]),rm_phix[j],
                                c(truncL[2*(k-1)+1],truncL[2*(k-1)+2]),truncQ[l]))
      }
    }
  }
}

#write.table( param, paste(pathRES, "/param.csv",sep=""), sep="\t" )
param<-as.data.frame(param)
colnames(param)<-c("maxEE","rmPhix","tcL","tcQ")

#######################################

irun = 1
table_taxa<-NULL
for (i in 1:nrow(param)){

	cat("**** Run ", irun, "\n")
	irun = irun + 1

  ####pipeline dada ########################################
  filtFs <- file.path(pathRES, "filtered2", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(pathRES, "filtered2", paste0(sample.names, "_R_filt.fastq.gz"))

  out<-filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=param$tcL[[i]],maxN=0, maxEE=param$maxEE[[i]], 
                     truncQ=param$tcQ[[i]], rm.phix=param$rmPhix[[i]],compress=TRUE, multithread = T,trimLeft =0)
  errF <- learnErrors(filtFs, multithread=TRUE)
  errR <- learnErrors(filtRs, multithread=TRUE)
  plotErrors(errF, nominalQ=TRUE)
  dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
  dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,trimOverhang = F, verbose=TRUE)
  seqtab <- makeSequenceTable(mergers)
  print(dim(seqtab))
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  taxa<- assignTaxonomy(seqtab.nochim, paste(pathDB,"/silva_nr_v132_train_set.fa.gz", sep=""), multithread=T,minBoot =50,outputBootstraps = F, tryRC=F )
  #taxa <- addSpecies(taxa, paste(pathDB, "silva_species_assignment_v132.fa.gz", sep="") )
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
  for(j in 1:nrow(taxa_data)){
    for(k in 2:7){
      nom_col<-paste(nom_col,taxa_data[j,k],sep="")
    }
    chemin_tax <-c(chemin_tax,nom_col)
    nom_col<-c()
  }
  ##ajout des colonnes chemin tax et du nombre de reads par asv 
  taxa_data<-cbind(taxa_data,chemin_tax,seqtab.nochim[1,])
  colnames(taxa_data)[ncol(taxa_data)]<-paste("track",i,sep="")
  colnames(taxa_data)[(ncol(taxa_data)-1)]<-c("nom")
  if(i>1){
    ##taxa
    
    #fusion des tables track i-1 et track i 
    table_fusion<-bind_rows(table_taxa,taxa_data)
    ##table comptabilisé les reads pour un chemin taxonomique
    table_count<-NULL  
    ##table taxonomique à construire
    table_tx<-NULL  
    name<-unique(table_fusion$nom)
    vec_sum<-c()
    for(j in 1:length(name)){
      ## recuperation de toutes les lignes ayant le meme chemin taxonomique
      table_count<-rbind(table_count,filter(table_fusion,nom==name[j]))
      ## on recupere les colonnes de la table taxonomique
      colonne_ta<-table_count[1,1:7]
      ## on somme les colonnes de reads de chaque track
      for(k in 8:ncol(table_fusion)){
        colonne_ta<-cbind(colonne_ta,sum(table_count[,k],na.rm=T))
      }
      ##on ajoute la nouvelle ligne du chemin taxonomique avec les reads sommés
      table_tx<-rbind(table_tx,colonne_ta)
      ## on réinitialise l'objet table_count 
      table_count<-NULL
    }
    table_taxa<-table_tx
    ##nom de colonnes
    for(k in 8:ncol(table_fusion)){
      colnames(table_taxa)[k]<-paste("track",(k-7),sep="")
    }
  }else{ 
    ##chemin taxonomique
    table_taxa<-taxa_data
  }
  rownames(out)<-NULL
  na<-nrow(taxa_na)
  if(is.null(nrow(taxa_na))){
    na<-0
  }
  track <- cbind(paste("track",i,sep=""),out,sum(getUniques(dadaFs)),sum(getUniques(dadaRs)),sum(getUniques(mergers)), rowSums(seqtab.nochim),na,length(unique(c(ind_e,ind))),
                 paste("Truncq_",param$tcQ[[i]],sep=""),paste("rmphix_",param$rmPhix[[i]],sep=""),paste("MaxEE_",paste(param$maxEE[[i]][1],param$maxEE[[i]][2],sep="-"),sep=""),
                 paste("TruncLen_",paste(param$tcL[[i]][1],param$tcL[[i]][2],sep="-"),sep=""))
  
  colnames(track) <- c("","input","filtered", "denoisedF", "denoisedR", "merged", "nonchim","assignes","badassignes","PRM_TRUNCQ","PRM_RMPHIX","PRM_MAXEE","Trunclen")
  ##changement du nom des colonnes de 
  for (j in 1:ncol(nochim)){
    colnames(nochim)[j]<-c(paste("asv",j,sep=""))
  }
  ##ecriture dans des fichiers 

  if(i>1){
    write.table(track,paste(sample.names,"_PDS.csv",sep=""),append=T,sep=",",qmethod="double",row.names = F,col.names = F)
  }else{
    write.table(track,paste(sample.names,"_PDS.csv",sep=""),append=F,sep=",",qmethod="double",row.names = F,col.names = T)
  }
}

asv<-c()  
for(j in 1:nrow(table_taxa)){
  asv<-c(asv,paste("asv",j,sep=""))
}
table_taxa<-cbind(asv,table_taxa)
colnames(table_taxa)[1]<-c("asv")
write.table(table_taxa[,1:7],paste(sample.names,"_taxa.csv",sep=""),append=F,sep=",",qmethod="double",row.names = F,col.names = T)
t<-cbind(table_taxa[,1],table_taxa[,9:ncol(table_taxa)])
write.table(t,paste(sample.names,"_counts.csv",sep=""),append=F,sep=",",qmethod="double",row.names = F,col.names = T)
T2<-Sys.time()
Tdiff<-difftime(T2, T1)

cat("End of pipeline\n")
