# Function that performs the first quality check: RLE
RLE_control<-function(Pset,n_probe_set,n_GSM,gsm,condizioni){

  # To identify GSMs that deviate from others using RLE, I can test the values returned by RLE using ANOVA.
  dati_RLE1<-RLE(Pset,type="value")
  n_probe_set<-nrow(dati_RLE1)
  dati_RLE=matrix(,nrow = n_probe_set*n_GSM,ncol=2)
  nomi=gsm
  dati_RLE[,2]=rep(nomi,each=n_probe_set)
  dati_RLE[,1]=as.vector(dati_RLE1)
  dati_RLE<-data.frame(espressione=dati_RLE[,1],gruppo=dati_RLE[,2])
  dati_RLE[,1]<-as.numeric(as.character(dati_RLE[,1]))
  
  # ANOVA
  res.aov_RLE<-aov(espressione~gruppo,data = dati_RLE)
  summ<-summary(res.aov_RLE)
  
  # If I have a p-value < 0.05, I reject the null hypothesis of equal means and proceed to identify which samples stand out.
  if (summ[[1]]$`Pr(>F)`[1]<0.05){
    # To check if there is any GSM that deviates significantly from all others and might need to be removed, I use the Tukey test.
    tuk<-TukeyHSD(res.aov_RLE)
    nomi_confronti<-rownames(tuk$gruppo)
    for (i in 1:n_GSM){
      indici=grep(nomi[i],nomi_confronti)
      p_val=tuk$gruppo[indici,4]
      # If a sample is found to deviate from all others in terms of mean, I classify it as a candidate for removal.
      if (sum(p_val<0.05)==length(indici))
      {
        condizioni[1,i]=1
      }else{
        condizioni[1,i]=0
      }
    }
  }else{
    condizioni[1,]=0
  }
  return(condizioni)
}

# Function that performs the second quality check: NUSE
NUSE_control<-function(Pset,condizioni,n_GSM){
  val=NUSE(Pset,type="value")
  soglia_NUSE=1.05  #from literature
  for (i in 1:n_GSM){
    N_med=median(as.numeric(val[,i]),na.rm = "TRUE")
    if(abs(N_med)>soglia_NUSE){
      condizioni[2,i]=1
    }else{
      condizioni[2,i]=0
    }
  }
  return(condizioni)
}

# Function that performs the third quality check on the distributions
distribution_control<-function(n_probes,n_GSM,nomi_gruppi,pmexp,condizioni){
  data_nuova=matrix(nrow = n_probes*n_GSM,ncol=2)
  data_nuova[,2]=rep(nomi_gruppi,rep(n_probes,each=n_GSM))
  data_nuova[,1]=as.vector(pmexp)
  data_nuova1<-data.frame(espressione=data_nuova[,1],gruppo=data_nuova[,2])
  # TEST KRUSKAL WALLIS
  data_nuova1$espressione<-as.numeric(as.character(data_nuova1$espressione))
  res.kw<-kruskal.test(formula=espressione~gruppo,data=data_nuova1)
  indici=vector()
  p_val=vector()
  
  # If p-value < 0.05, the null hypothesis of equal means is rejected, and the 
  # different GSMs are examined.
  if (res.kw$p.value<0.05){
    # Tukey test to check which samples consistently deviate (2vs2)
    tuk<-pairwise.wilcox.test(data_nuova1$espressione,data_nuova1$gruppo,p.adjust.method = "bonferroni")
    nomi_confronti<-colnames(tuk$p.value)
    p_val=vector()
    for (i in 1:n_GSM){
      if (i==1){
        p_val=tuk$p.value[,i]
      }else if (i==n_GSM){
        p_val=tuk$p.value[i-1,]
      }else{
        p_val=c(tuk$p.value[,i],tuk$p.value[i-1,])
      }
      p_val=na.omit(p_val)
      # If the sample is found to deviate from all others in terms of mean, I classify it as a candidate for removal.      
      if (sum(p_val<0.05)==length(nomi_gruppi)-1){
        condizioni[3,i]=1
      }else{
        condizioni[3,i]=0
      }
    }
  }else{
    condizioni[3,]=0
  }
  return(condizioni)
}

# Function that identifies the GSMs to be eliminated
GSM_delete<-function(condizioni){
  da_eliminare=vector()
  for (i in 1:length(condizioni[1,])){
    # If an array fails to meet any of the three quality checks, it will be eliminated.
    if (sum(condizioni[,i])==nrow(condizioni)){
      da_eliminare[i]=1
    }else{
      da_eliminare[i]=0
    }
  }
  da_eliminare<-which(da_eliminare>0)
  return(da_eliminare)
}

# Function that removes low quality sample
post_eliminazione<-function(da_eliminare,nomi_gruppi,data){

  if (length(da_eliminare)!=0){
    nomi_gruppi<-nomi_gruppi[-da_eliminare]
    data<-data[,-da_eliminare]
  }
  nomi_gruppi_unici<-unique(nomi_gruppi)
  n_gruppi<-length(nomi_gruppi_unici)
  n_repliche<-c()
  for(i in 1:n_gruppi){
    n_repliche[i]<-length(which(nomi_gruppi_unici[i]==nomi_gruppi))
  }
  
  risultati<-matrix(,nrow = 2,ncol = length(nomi_gruppi_unici))
  risultati[1,]<-n_repliche
  risultati[2,]<-nomi_gruppi_unici
  
  output<-list(risultati,data)
  return(output)
}


# Function that analyzes .CEL files --> Affymetrix
affymetrix_analysis<-function(tmp_MA_path,GSE,k,locus_tag,locus_tag_unici,gsm,n_GSM,group_GSM,metodo_correzione,PVALUE){
  
  celpath = paste(tmp_MA_path,GSE[k],sep="/")
  celFiles <- list.celfiles(celpath, full.names=TRUE)
  # Read .CEL files
  data<-read.celfiles(celFiles)

  pvalue<-c()
  # Number of probes
  n_probes<-length(probeNames(data))
  # Number of probeset
  n_probe_set<-length(featureNames(data))

  #### 1.2) PSEUDO RECONSTRUCTION OF THE CHIP
  Pset=fitProbeLevelModel(data)
  Pset=na.omit(Pset)

  # A matrix is initialized that will contain either 1 or 0 depending on whether 
  # the conditions are met or not. (1 indicates the condition is not met, 0 
  # indicates the condition is met).
  # By summing the columns, it can be determined how many samples do not meet 
  # the conditions and are candidates for elimination.  
  condizioni=matrix(,nrow=3,ncol=n_GSM)
  
  #### 1.3) RLE (Relative Logarithmic Expression)
  condizioni<-RLE_control(Pset,n_probe_set,n_GSM,gsm,condizioni)

  #### 1.4) NUSE (Normalized Unscaled Standard Error)
  condizioni<-NUSE_control(Pset,condizioni,n_GSM)
  
  #### 2) DISTRIBUTIONS OF EXPRESSION VALUES
  gruppi_unici<-unique(group_GSM)
  n_gruppi<-length(gruppi_unici)
  n_repliche<-c()
  for(i in 1:length(gruppi_unici)){
    n_repliche[i]<-length(which(gruppi_unici[i]==group_GSM))
  }
  # perfect match expression
  pmexp=pm(data)

  condizioni<-distribution_control(n_probes,n_GSM,gsm,pmexp,condizioni)
  #### GSM TO BE REMOVED
  da_eliminare<-GSM_delete(condizioni)
  # Number of replicates for each group after the removal of samples
  output<-post_eliminazione(da_eliminare,group_GSM,data)
  risultati<-output[[1]]
  n_repliche<-risultati[1,]
  nomi_gruppi_unici<-risultati[2,]
  new_data<-output[[2]]
  
  #### 3) NORMALIZATION with RMA (Robust Multichip Average)
  data.rma=rma(new_data)
  
  n_GSM_post_eliminazione<-length(new_data@phenoData@data$index)
  gsm_post_elim<-gsm[new_data@phenoData@data$index]

  # 4) DISTINCTION BETWEEN DE AND NON-DE GENES
  # The expression of each gene is derived from the expression of all probesets 
  # as the average of the probesets related to that gene.
  espressione_gene<-matrix(nrow = length(locus_tag_unici),ncol = n_GSM_post_eliminazione)
  expr<-exprs(data.rma)
  for (i in 1:length(locus_tag_unici)){
    ind<-which(locus_tag_unici[i]==locus_tag)
    if(length(ind)>1){
      for (j in 1:n_GSM_post_eliminazione){
        espressione_gene[i,j]<-mean(expr[ind,j])
      }
    }else{
      espressione_gene[i,]<-expr[i,]
    }
  }
  
  #### 4.1 KRUSKAL WALLIS TRA GRUPPI
  # An ANOVA test needs to be performed for each gene to assess whether its expression varies based on group membership.
  data_nuova_post_rma=matrix(nrow = n_GSM_post_eliminazione,ncol = 2)
  nomi_gruppi_post_rma<-nomi_gruppi_unici
  data_nuova_post_rma[,2]=rep(unique(nomi_gruppi_post_rma),n_repliche)
  
  # Kruskal-Wallis test
  for (t in 1:length(locus_tag_unici)){
    data_nuova_post_rma[,1]=as.vector(espressione_gene[t,])
    data_nuova_post_rma1<-data.frame(espressione=data_nuova_post_rma[,1],gruppo=data_nuova_post_rma[,2])
    data_nuova_post_rma1$espressione<-as.numeric(as.character(data_nuova_post_rma1$espressione))
    res.kw_post_rma<-kruskal.test(formula=espressione~gruppo,data=data_nuova_post_rma1)
    pvalue[t]<-res.kw_post_rma$p.value
    if(pvalue[t]=="NaN"){
      pvalue[t]<-0
    }
  }
  
  # pvalue correction
  if (metodo_correzione %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) {
  	p_adjust <- p.adjust(pvalue, method = metodo_correzione, n = length(locus_tag_unici))
	} else {
  	warning("The input correction method is not valid: choose between: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' or 'none'.\nUsing 'none' by default..")
		p_adjust <- p.adjust(pvalue, method = "none",  n = length(locus_tag_unici))
	}

  ind<-which(p_adjust>PVALUE)
  NDE<-rep(0,length(p_adjust))
  NDE[ind]<-1

  results<-cbind(NDE,locus_tag_unici,espressione_gene)
  colnames(results)<-c("NON_DE","locus_tag_unici",gsm_post_elim)
  
  return(results)
}

