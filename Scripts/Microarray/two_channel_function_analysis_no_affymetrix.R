# Function that reads two-channel files 
# INPUT: tipo_scan, the list of files to analyze, the path where the files are located, 
# the number of channels
# OUTPUT: RG --> structure containing the expression data of the two channels
lettura_file_2canali<-function(tipo_scan,elenco_file,GSM_path){
  RG<-c()
  setwd(GSM_path)
  stringsFG<-c('[Rr]aw','F[0-9]+.[mM]ean$','[rg][mM]ean[sS]ignal','[rg][mM]edian[sS]ignal','Signal [mM]ean_?[cC]y[0-9]')
  stringsBG<-c('[Bb]ackground','B[0-9]+.[mM]ean$','[rg]BG[mM]ean[sS]ignal','[rg]BG[mM]e[sS]ignal','Background [mM]ean_?[cC]y[0-9]')
  lines<-readLines(elenco_file[1])
  split_lines<-strsplit(lines,"\t")
  for(i in 1:length(split_lines)){
    indexBG<-c()
    indexFG<-c()
    for(j in 1:length(stringsFG)){
      indexFG<-which(regexpr(stringsFG[j],gsub("\"","",split_lines[[i]]))>0)
      indexBG<-which(regexpr(stringsBG[j],gsub("\"","",split_lines[[i]]))>0)
      if(length(indexFG)>0 && length(indexBG)>0){
        verde<-gsub("\"","",split_lines[[i]][indexFG[1]])
        rosso<-gsub("\"","",split_lines[[i]][indexFG[2]])
        verdeBG<-gsub("\"","",split_lines[[i]][indexBG[1]])
        rossoBG<-gsub("\"","",split_lines[[i]][indexBG[2]])
        break
      }
    }
    if(length(indexFG)>0 || length(indexBG)>0){
      break
    }
  }
 
  for (i in 1:length(tipo_scan)){
      out<-tryCatch(
        expr={
          RG<-read.maimages(files = elenco_file, source = tipo_scan[i],path = GSM_path,columns = list(G=verde,R=rosso,Gb=verdeBG,Rb=rossoBG))
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
        }
      )
    
    if (length(RG)>0){
      break;
    }
  }
  
  return(RG)
}

# M/A ratio
# INPUT: the number of GSMs, the vector containing the unique locus tags, 
# the vector with all the locus tags, the structure containing the expression 
# data after between-array normalization
# OUTPUT: the matrix containing the M/A ratio and A_value
calcolo_rapporto_M_su_A<-function(n_GSM,locus_tag_unici,locus,MA.pAq){
  rapporto_M_su_A<-matrix(ncol=n_GSM,nrow = length(locus_tag_unici))
  A_value<-matrix(ncol=n_GSM,nrow = length(locus_tag_unici))
  for (t in 1:n_GSM){
    for (r in 1:length(locus_tag_unici)){
      ind<-which(locus==locus_tag_unici[r])
      M_value<-mean(MA.pAq$M[ind,t])
      A_value[r,t]<-mean(MA.pAq$A[ind,t])
      rapporto_M_su_A[r,t]<-M_value/A_value[r,t]
    }
  }
  output<-list(rapporto_M_su_A,A_value)
  return(output)
}

# Function that reads two-channel files
# INPUT: scanner type, list of files to be analyzed, path where the files are 
# located, number of GSMs, vector of unique locus tags, vector with all locus
# tags, vector with GSM codes, number of GSE, # index of the analyzed GSE, statistical test correction method, 
# GSE codes, GSM groups
# OUTPUT: matrix where the first column represents DE/NON DE genes and the 
# subsequent columns are the expression values
two_channel_analysis<-function(tipo_scan,elenco_file,GSM_path,n_GSM,locus_tag_unici,locus,gsm,metodo_correzione, PVALUE){
  # Leggo i file e li salvo in RG: questa struttura conterr? i dati di verde e rosso dei due canali per ogni array
  RG<-lettura_file_2canali(tipo_scan,elenco_file,GSM_path)
  
  if(is.null(RG)==0){
    
    #### QUALITY CONTROL AND NORMALIZATION
    RG.b <- limma::backgroundCorrect(RG, method="normexp", offset=1)
    MA.p <- normalizeWithinArrays(RG.b, method="loess")
    MA.pAq <- normalizeBetweenArrays(MA.p, method="Aquantile")

    # Calcolo il rapporto M/A
    output<-calcolo_rapporto_M_su_A(n_GSM,locus_tag_unici,locus,MA.pAq)
    rapporto_M_su_A<-output[[1]]
    A_value<-output[[2]]
    
    # TEST DI WILCOXON
    pvalue<-c()
    for (i in 1:length(locus_tag_unici)){
      rr<-wilcox.test(rapporto_M_su_A[i,],conf.level = 0.05)
      pvalue[i]<-rr$p.value
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
    
    results<-cbind(NDE,locus_tag_unici,A_value)
    colnames(results)<-c("NON_DE","locus_tag_unici",gsm)
  }else{
    results<-c()
  }
    return(results)
}

