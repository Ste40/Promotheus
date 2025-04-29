# Function to read files
# INPUT: scanner type used, list of files, path where the files are located
# OUTPUT: green channel expression data for each GSM
lettura_file_1canale<-function(tipo_scan,elenco_file,GSM_path){
  G<-c()
  setwd(GSM_path)
  stringsFG<-c('[Rr]aw','F[0-9]+.[mM]ean$','[g][mM]ean[sS]ignal','[g][mM]edian[sS]ignal','^PM$')
  stringsBG<-c('[Bb]ackground','B[0-9]+.[mM]ean$','[g]BG[mM]ean[sS]ignal','[g]BG[mM]e[sS]ignal','^MM$')
  lines<-readLines(elenco_file[1])
  split_lines<-strsplit(lines,"\t")
  for(i in 1:length(split_lines)){
    ind<-c()
    for(j in 1:length(stringsFG)){
      indexFG<-which(regexpr(stringsFG[j],gsub("\"","",split_lines[[i]]))>0)
      indexBG<-which(regexpr(stringsBG[j],gsub("\"","",split_lines[[i]]))>0)
      if(length(indexFG)>0){
        verde<-gsub("\"","",split_lines[[i]][indexFG[1]])
      }
      if(length(indexBG)>0){
        verdeBG<-gsub("\"","",split_lines[[i]][indexBG[1]])
        break
      }
    }
  }
  
  for (i in 1:length(tipo_scan)){
      out<-tryCatch(
        expr={
          G<-read.maimages(files = elenco_file, source = tipo_scan[i],path = GSM_path,green.only = T,columns = list(G=verde,Gb=verdeBG))
        },
        error=function(cond) {
          message(paste("Errore nella lettura dei file", tipo_scan[i],"\n"))
          message("Messaggio di errore originario:")
          message(paste(cond,"\n"))
          return(NA)
        },
        warning=function(cond) {
          message(paste("Warning nella lettura dei file:", tipo_scan[i],"\n"))
          message("Messaggio di warning originario:")
          message(paste(cond,"\n"))
          return(NULL)
        },
        finally={
          message(paste("Tipo scan processato:", tipo_scan[i],"\n"))
          message("Some other message at the end")
        }
      )
      if (length(G)>0){
        break;
      }
    
  }
  return(G)
}

# Funzione che effettua un mapping delle sonde (potrebebro esserci sonde di controllo non relative a geni)
# INPUT: path in cui si trovano i file dei GSM, elenco dei file, sonde relative alla GPL corrente
# OUTPUT: indici relativi alle sonde da considerare  -----> CAMBIATO
mapping_sonde<-function(GSM_path,elenco_file,sonde_GPL,vet_sonde){
  # Leggo il testo ed elimino le righe da ignorare
  testo<-read.table(paste(GSM_path,elenco_file[1],sep="/"),sep="\t",blank.lines.skip = FALSE,quote="") # aggiunto quote=""
  ind<-max(which(testo[,2]==""))
  ind1<-max(which(testo[,1]=="TYPE"))
  if (ind>ind1 & ind!=length(testo[,1])){
    colnames(testo)<-testo[ind+1,]
    for (j in 1:(ind+1)){
      testo<-testo[-1,]
    }
  }else if(ind1>ind){
    colnames(testo)<-testo[ind1+1,]
    for (j in 1:(ind1+1)){
      testo<-testo[-1,]
    }
  }else{
    colnames(testo)<-testo[1,]
    testo<-testo[-1,]
  }
  # Recupero la colonna relativa alle sonde riportate nel GSM
  coln<-colnames(testo)
  
  ind<-which(coln=="ProbeName")
  if(length(ind)==0){
    ind<-which(coln=="SEQ_ID")
  }
  if(length(ind)==0){
    ind<-which(coln=="PROBE_ID")
  }
  sonde_GSM<-testo[,ind]
  if(is.na(vet_sonde[1])==0){
    sonde_GSM<-vet_sonde
  }
  ind<-which(is.na(sonde_GPL))
  if(length(ind)>0){
    sonde_GPL<-sonde_GPL[-ind]
  }
  # Ottengo gli indici delle sonde da tenere
  ind1<-matrix(,ncol = length(sonde_GPL),nrow=100)
  for(i in 1:length(sonde_GPL)){
    ind2<-which(sonde_GPL[i]==sonde_GSM)
    if(length(ind2)>0){
      ind1[1:length(ind2),i]<-ind2
    }
    
  }
  ind<-as.vector(ind1)
  ind<-ind[-which(is.na(ind)==1)]
  ind<-unique(ind)
  return(ind)
}

# Function to analyze one-channel chip
one_channel_analysis<-function(n_GSM, tipo_scan,elenco_file,GSM_path,locus_tag_unici,sonde_GPL,locus_tag,metodo_correzione,group_GSM,vet_sonde,PVALUE,gsm){
  
  G<-lettura_file_1canale(tipo_scan,elenco_file,GSM_path)
  if(length(G)>0){
    # BACKGROUND CORRECTION AND NORMALIZATION
    Gb <- limma::backgroundCorrect(G, method="normexp", offset=1)
    Gn <- normalizeBetweenArrays(Gb, method="quantile")

    # Mapping probes
    ind<-mapping_sonde(GSM_path,elenco_file,sonde_GPL,vet_sonde)
    if (length(ind)>1){
      esp<-Gn$E[ind,]
      vet_sonde<-vet_sonde[ind]
    }else{
      esp<-Gn$E
    }
    # Expression values
    espressione_geni<-matrix(,ncol=n_GSM,nrow = length(locus_tag_unici))
    for (i in 1:length(locus_tag_unici)){
      ind<-which(regexpr(locus_tag_unici[i],vet_sonde)==1)
      if (length(ind)==1){
        espressione_geni[i,]<-as.numeric(esp[ind,])
      }else if (length(ind)>1){
        for (j in 1:n_GSM){
          espressione_geni[i,j]<-mean(as.numeric(esp[ind,j]))
        }
      }
    }
    
    pvalue<-c()
    nomi_gruppi_ordinati<-group_GSM
    data<-matrix(ncol = 2,nrow = n_GSM)
    data[,2]<-nomi_gruppi_ordinati
    # Kruskal-Wallis test
    for (t in 1:length(locus_tag_unici)){
      data[,1]=as.vector(espressione_geni[t,])
      data_DF<-data.frame(espressione=data[,1],gruppo=data[,2])
      data_DF$espressione<-as.numeric(as.character(data_DF$espressione))
      res.kw_post_rma<-kruskal.test(formula=espressione~gruppo,data=data_DF)
      pvalue[t]<-res.kw_post_rma$p.value
      if(is.na(pvalue[t])==1){
        pvalue[t]<-1
      }
    }
    
    if (metodo_correzione %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) {
      p_adjust <- p.adjust(pvalue, method = metodo_correzione, n = length(locus_tag_unici))
    } else {
      warning("The input correction method is not valid: choose between: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' or 'none'.\nUsing 'none' by default..")
      p_adjust <- p.adjust(pvalue, method = "none",  n = length(locus_tag_unici))
    }
    ind<-which(p_adjust>PVALUE)
    NDE<-rep(0,length(p_adjust))
    NDE[ind]<-1
    
    results<-cbind(NDE,locus_tag_unici,espressione_geni)
    colnames(results)<-c("NON_DE","locus_tag_unici",gsm)
  }else{
    results<-c()
  }
  
  return(results)
}
