# Function that returns the file extension for the current GPL
control_extension_file<-function(elenco_file){
  expression<-c()
  expression[1]<-".txt"
  expression[2]<-".CEL"
  expression[3]<-".pair"
  expression[4]<-".gpr"
  ind<-c()
  for(i in 1:length(expression)){
    ind[i]<-regexpr(expression[i],elenco_file[1])>0
  }
  ind<-which(ind==1)
  return(expression[ind])
}

# Function that maps the genes used in the analyzed experiment to those in the reference genome
mapping<-function(risultati, locus_tag_unici, geni_ortologhi, gsm, path_file, GSE_code){
  
  for(i in 1:length(locus_tag_unici)){
    ind_col<-which(locus_tag_unici[i]==geni_ortologhi,arr.ind=T)
    if(length(ind_col)==2){
      ind_col<-as.numeric(ind_col[2])
      break
    }else if(length(ind_col)>2){
      ind_col<-as.numeric(ind_col[1,2])
      break
    }
  }
  riga<-c()
  for(i in 1:length(locus_tag_unici)){
    ind<-which(locus_tag_unici[i]==geni_ortologhi[,ind_col],arr.ind=T)
    if(length(ind)>0){
      riga[i]<-ind[1]
    }else{
      riga[i]<-NA
    }
  }
  ind<-which(is.na(riga))
  if(length(ind)>0){
    riga<-riga[-ind]  
    risultati1<-risultati[-ind,]
    locus_tag_unici<-locus_tag_unici[-ind]
  }
  
  DE<-risultati1[,1]
  mean_ex<-c()
  sd_ex<-c()
  cv_ex<-c()
  gene_name<-geni_ortologhi[riga,1]
  locus_tag_unici<-geni_ortologhi[riga,2]
  espr<-risultati1[,3:ncol(risultati1)]
  for(i in 1:length(riga)){
    mean_ex[i]<-mean(as.numeric(espr[i,]))
    sd_ex[i]<-sd(as.numeric(espr[i,]))
    cv_ex[i]<-sd_ex[i]/mean_ex[i]
  }
  
  # percentiles
  percentile<-matrix(ncol = ncol(espr),nrow=length(locus_tag_unici))
  for(i in 1:ncol(percentile)){
    rankin<-rank(as.numeric(espr[,i]))
    percentile[,i]<-rankin/length(locus_tag_unici)*100
  }
  # mean percentiles
  percentile_medio<-c()
  for(i in 1:nrow(percentile)){
    percentile_medio[i]<-mean(percentile[i,])
  }
  
  final_table<-cbind(gene_name,locus_tag_unici,espr,mean_ex,sd_ex,cv_ex,DE,percentile,percentile_medio)
  colnames(final_table)<-c("Geneid","locus_tag",gsm, "means","sd","cv","diff_signif",c(1:length(gsm)), "Percentiles_mean")
  final_table<-data.frame(final_table)
  fwrite(final_table,paste(path_file,paste(GSE_code,"_full_tab.txt",sep=""),sep="/"),sep="\t",row.names = F)
}

# Function to analyze each experiment
intra_exp_analysis<-function(tmp_MA_path, tmp_genome_path, CORRECTION_DE, PVALUE){
  
  tab1<-read.xlsx(file.path(tmp_MA_path,"total_metadata_GSE.xlsx"))
  GSE<-tab1$GSE
  GPL<-tab1$GPL
  organism<-tab1$organism
  taxid<-tab1$taxid
  n_canali<-tab1$n_canali
  n_GSM<-tab1$n_GSM
  manufacturer<-tab1$manufacturer
  asm_name<-tab1$asm_name
  n_GSE<-length(GSE)
  GPL_uniche<-unique(GPL)
  taxid_species<-tab1$taxid_species
  
  GSM<-read.xlsx(file.path(tmp_MA_path,"total_GSM.xlsx"))
  
  gsm_condition<-read.xlsx(file.path(tmp_MA_path,"gsm_condition.xlsx"))
  
  tabella_mapping_sonde<-read.xlsx(file.path(tmp_MA_path,"total_probes.xlsx"))
  
  geni_ortologhi<-read.table(file.path("/output/orthologous_genes.txt"))

  # types of scanners that can be used
  tipo_scan<-c("generic","agilent","arrayvision","bluefuse","genepix","imagene","quantarray","scanarrayexpress","smd","spot")
  
  for (k in 1:n_GSE){  
    mex<-paste("######### GSE ANALYSED:", GSE[k],"#########",sep=" ")
    message(mex)
    
    GPL_file<-paste(GPL[k],".txt",sep="")
    data<-read.table(file.path(tmp_MA_path,GPL_file),header = T,sep="\t",fill = T)
    sonde_GPL<-data$ID_probe
    locus<-data$locus
    gene_name<-data$gene_name
    
    locus_tag_unici<-unique(locus)
    ind<-which(is.na(GSM[,k]))
    if(length(ind)>0){
      gsm<-GSM[-ind,k]
    }else{
      gsm<-GSM[,k]
    }
    
    n_GSM<-length(gsm)
    group_GSM<-gsm_condition$group[which(gsm_condition$GSE==GSE[k])]
    group_GSM<-strsplit(group_GSM,",")[[1]]
    
    path<-paste(tmp_MA_path,GSE[k],sep="/")
    elenco_file<-dir(path,pattern = "^GSM")
    extension<-control_extension_file(elenco_file)
    if(length(extension)!=0){
      if(regexpr(".CEL",extension)==1){
        risultati<-affymetrix_analysis(tmp_MA_path,GSE,k,locus,locus_tag_unici,gsm,n_GSM,group_GSM,CORRECTION_DE, PVALUE)
      }else if(n_canali[k]==2){
        risultati<-two_channel_analysis(tipo_scan,elenco_file,path,n_GSM,locus_tag_unici,locus,gsm,CORRECTION_DE,PVALUE)
      }else{
        risultati<-one_channel_analysis(n_GSM,tipo_scan,elenco_file,path,locus_tag_unici,sonde_GPL,locus,CORRECTION_DE,group_GSM,tabella_mapping_sonde[,k], PVALUE,gsm)
      }
    }

    
    if(length(extension)==0 || length(risultati)==0){
      mex<-paste("######### Error during the analysis of", GSE[k],"- Files not read #########",sep=" ")
      message(mex)
      testo0<-paste('bash -c',shQuote(paste("rm -r ",tmp_MA_path,"/",GSE[k],sep="")),sep=" ")
      system(testo0)
    }else{
      locus_tag_unici<-risultati[,2]
      gsm_post_elim<-colnames(risultati)[3:ncol(risultati)]
      mapping(risultati, locus_tag_unici, geni_ortologhi, gsm_post_elim, tmp_MA_path, GSE[k])

    }
  }
  
  GSE_finals<-dir(tmp_MA_path,"full_tab")
  if(length(GSE_finals)==0){
    message(" ################## It was not possible to analyze any of the selected experiments. ########################")
  }
}
