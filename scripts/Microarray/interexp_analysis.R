################################################################################
########################### INTER EXPERIMENT ANALYSIS ##########################
################################################################################
# Merge microarray and RNAseq files
unique_matrix<-function(tmp_MA_path, path_rna,tmp_genome_path){
  
  geni_ortologhi<-read.table("/output/orthologous_genes.txt")
  #MA
  file_MA<-dir(tmp_MA_path,pattern="_full_tab")
  gene_com<-c()
  tab<-list()
  o<-1
  GSEcode<-c()
  if(length(file_MA)>0){
    for(i in 1:length(file_MA)){
      tab[[i]]<-read.table(file.path(tmp_MA_path,file_MA[i]),header = T,sep = "\t")
      if(i==1){
        gene_com<-tab[[i]]$locus_tag
      }else{
        gene_com<-intersect(tab[[i]]$locus_tag,gene_com) 
      }
    }
    o<-length(file_MA)+1
    GSEcode<-gsub("_full_tab.txt","",file_MA)
  }

  #RNAseq
  file_RNAseq<-dir(path_rna,pattern="full_tab")
  if(length(file_RNAseq)>0){
    for(i in 1:length(file_RNAseq)){
      tab[[o]]<-read.table(file.path(path_rna,file_RNAseq[i]),header = T,sep="\t")
      locus_tag_unici<-tab[[o]]$locus_tag
      
      for(x in 1:length(locus_tag_unici)){
        ind_col<-which(locus_tag_unici[x]==geni_ortologhi,arr.ind=T)
        if(length(ind_col)==2){
          ind_col<-as.numeric(ind_col[2])
          break
        }else if(length(ind_col)>2){
          ind_col<-as.numeric(ind_col[1,2])
          break
        }
      }
      riga<-c()
      for(x in 1:length(locus_tag_unici)){
        ind<-which(locus_tag_unici[x]==geni_ortologhi[,ind_col],arr.ind=T)
        if(length(ind)>0){
          riga[x]<-ind[1]
        }else{
          riga[x]<-NA
        }
      }
      ind<-which(is.na(riga))
      if(length(ind)>0){
        riga<-riga[-ind]  
        tab[[o]]<-tab[[o]][-ind,]
        locus_tag_unici<-locus_tag_unici[-ind]
        tab[[o]]$locus_tag<-geni_ortologhi[riga,2]
      }

      gene_com<-intersect(tab[[o]]$locus_tag,gene_com) 
      o<-o+1
    }
    GSEcode<-c(GSEcode,gsub("full_tab.txt","",file_RNAseq))
  }
  

  if(OPERONS_DETECTION){
    directories<-dir("/output/RNA_Seq")
    ind<-which(directories=="results")
    if(length(ind)>0){
      directories<-directories[-ind]
    }
    ind<-which(directories=="Experimental_Groups.xlsx")
    if(length(ind)>0){
      directories<-directories[-ind]
    }
    for(i in 1:length(directories)){
      GSE_fulltab<-dir("/output/RNA_Seq/results")
      ind<-which(regexpr(directories[i],GSE_fulltab)>0)
      if(length(ind)>1){
        ind<-ind[1]
      }
      GSE_fulltab<-read.table(file.path("/output/RNA_Seq/results",GSE_fulltab[ind]),sep="\t",header = TRUE)
      locus_fulltab<-GSE_fulltab$locus_tag
      tag<-GSE_fulltab$n.genes

      ind<-which(tag=="not first gene in operon")
      if(length(ind)>0){
        locus_fulltab<-locus_fulltab[-ind]
      }

      for(x in 1:length(locus_fulltab)){
        ind_col<-which(locus_fulltab[x]==geni_ortologhi,arr.ind=T)
        if(length(ind_col)==2){
          ind_col<-as.numeric(ind_col[2])
          break
        }else if(length(ind_col)>2){
          ind_col<-as.numeric(ind_col[1,2])
          break
        }
      }

      riga<-c()
      for(x in 1:length(locus_fulltab)){
        ind<-which(locus_fulltab[x]==geni_ortologhi[,ind_col],arr.ind=T)
        if(length(ind)>0){
          riga[x]<-ind[1]
        }else{
          riga[x]<-NA
        }
      }

      ind<-which(is.na(riga))
      if(length(ind)>0){
        riga<-riga[-ind]  
        locus_fulltab<-geni_ortologhi[riga,2]
      }
      if(length(locus_fulltab)>0){
        gene_com<-intersect(locus_fulltab,gene_com) 
        if(length(gene_com)==0){
          message("ERROR: all genes common to all experiments are organized in operons, therefore no results can be generated.")
        }
      }
    }
  }

  if(length(gene_com)>0){
    #Union expression
    ESP<-c()
    n_GSM<-c()
    DE<-rep("NO",length(gene_com))
    for(i in 1:length(tab)){
      ind<-c()
      for(j in 1:length(gene_com)){
        ind[j]<-which(gene_com[j]==tab[[i]]$locus_tag)
      }
      ind_mean<-which(colnames(tab[[i]])=="means")
      if(length(ESP)==0){
        ESP<-as.matrix(tab[[i]][ind,3:(ind_mean-1)])
      }else{
        ESP<-cbind(ESP,as.matrix(tab[[i]][ind,3:(ind_mean-1)]))
      }
      n_GSM[i]<-ind_mean-3
       # nolint: trailing_whitespace_linter.
      # DE
      inde_DE<-which(colnames(tab[[i]])=="diff_signif")
      ind_YES<-which(tab[[i]][ind,inde_DE]=="YES")
      if(length(ind_YES)>0){
        DE[ind_YES]<-"YES"
      }
    }

    gene_name<-c()
    for(i in 1:length(gene_com)) {
      ind <- which(gene_com[i] == geni_ortologhi,arr.ind=TRUE)
      if(length(ind)==2){
        gene_name[i] <- geni_ortologhi[ind[1],1]
      }else if(length(ind)>2){
        gene_name[i] <- geni_ortologhi[ind[1,1],1]
      }
      
    }
    
    batches<-c(rep(GSEcode,n_GSM))
    colnames(ESP)<-batches
    final_tab<-cbind(Geneid=gene_name,locus_comuni=gene_com,ESP,DE)
    output<-list(final_tab,batches)
  }else{
    output<-list()
  }

  return(output)
}

# Function to create final table
create_final_table<-function(espressione_media,tabella_completa){
  
  gene_name<-tabella_completa[,1]
  locus_tag<-tabella_completa[,2]
  DE<-tabella_completa[,ncol(tabella_completa)]
  
  elenco_GSE<-colnames(espressione_media)
  
  mean_psnorm<-c()
  sd_psnorm<-c()
  for(i in 1:length(locus_tag)){
    sd_psnorm[i]<-sd(espressione_media[i,])
    mean_psnorm[i]<-mean(espressione_media[i,])
  }
  
  # percentiles
  percentile<-matrix(ncol = ncol(espressione_media),nrow=nrow(espressione_media))
  for(i in 1:ncol(percentile)){
    rankin<-rank(as.numeric(espressione_media[,i]))
    percentile[,i]<-rankin/nrow(espressione_media)*100
  }
  percentile_medio<-c()
  percentile_sd<-c()
  for(i in 1:nrow(percentile)){
    percentile_medio[i]<-mean(percentile[i,])
    percentile_sd[i]<-sd(percentile[i,])
  }
  
  # final table
  finaltab<-cbind(gene_name,locus_tag,espressione_media,mean_psnorm,sd_psnorm,percentile,percentile_medio,percentile_sd,DE)
  etichette_espressione<-paste('Expression_',elenco_GSE,sep="")
  etichette_percentile<-paste('Percentile_',elenco_GSE,sep="")
  colnames(finaltab)<-c("Gene_name","GeneID",etichette_espressione,"Mean_Expression","CV",etichette_percentile,"Mean_Percentile","Percentile_SD","DE")
  finaltab<-data.frame(finaltab)
  
  return(finaltab)
}

# Remove temporary files 
last_operations<-function(tmp_MA_path,path_rna,path_fulltab_res,path_output){

  files <- list.files(tmp_MA_path, pattern = "full_tab", full.names = TRUE)
  path_dest <- path_fulltab_res
  destination_paths <- file.path(path_dest, basename(files))
  file.rename(files, destination_paths)
  
  files <- list.files(path_rna, pattern = "full_tab", full.names = TRUE)
  path_dest <- path_fulltab_res
  destination_paths <- file.path(path_dest, gsub(" ", "_", basename(files)))
  file.copy(files, destination_paths)

  setwd(path_output)
  testo0<-paste('bash -c',shQuote(paste("rm -r tmp",sep=" ")),sep=" ")
  system(testo0)
}

# Function to add promter sequece (300bp) to final table
add_promoter_sequence<-function(path_output,finaltab,asm_pref){
  summary<- fread("/output/assembly_summary.txt", quote=F)
  
  ind<-which(summary$asm_name==asm_pref)
  ftp_path<-summary$ftp_path[ind]
  ftp_string<- gsub(".*/", "", ftp_path)
  GFF_ftp<- paste0(ftp_path, "/", ftp_string, "_genomic.gff.gz")
  fna_ftp<-paste0(ftp_path, "/", ftp_string, "_genomic.fna.gz")
  if(file.exists(paste(path_output,"/",asm_pref,".gff",sep=""))==0){
    for(i in 1:3){
      tryCatch({
        annot<- curl_download(GFF_ftp, destfile=paste0(path_output,"/", ftp_string, "_genomic.gff.gz"), mode="wb")
        if(length(annot)>0){
          break
        }
      },error=function(e){
        message("An error occurred while downloading .gff file")
        Sys.sleep(60)
      })
    }
    if(length(annot)>0){
      new_name<-gsub(".+gff.gz",paste(path_output,"/",asm_pref,".gff",sep=""),annot)
      gunzip(filename=annot,destname=new_name,overwrite=T)
    }
    
  }
  
  if(file.exists(paste(path_output,"/",asm_pref,".fna",sep=""))==0){
     for(i in 1:3){
      tryCatch({
        annot<- curl_download(fna_ftp, destfile=paste0(path_output,"/", ftp_string, "_genomic.fna.gz"), mode="wb")
        if(length(annot)>0){
          break
        }
      },error=function(e){
        message("An error occurred while downloading .gff file")
        Sys.sleep(60)
      })
    }
    if(length(annot)>0){
      new_name<-gsub(".+fna.gz",paste(path_output,"/",asm_pref,".fna",sep=""),annot)
      gunzip(filename=annot,destname=new_name,overwrite=T)
    }
    
  }

  GFF_file<-file.path(paste(path_output, "/", asm_pref, ".gff", sep=""))
  fna_file<-file.path(paste(path_output, "/", asm_pref, ".fna", sep=""))
  if(length(GFF_file)>0 & length(fna_file)>0){
    GFF<-read.gff(file.path(GFF_file))
    FASTA<-read.fasta(file.path(fna_file))[[1]]
    #FASTA<-FASTA[[1]]
  
    N<-300
    Sequence<-c()
    for(i in 1:nrow(finaltab)){
      ind<-which(regexpr(finaltab$GeneID[i],GFF$V9)>0)
      ind<-ind[1]
      start<-GFF[ind,4]
      end<-GFF[ind,5]
      strand<-as.character(GFF[ind,7])
    
      if(strand=="-"){
        start_s<-end+1
        end_s<-start_s+N-1
      }else{
        end_s<-start-1
        start_s<-end_s-N+1
      }
    
      if(start_s>0 && end_s<length(FASTA)){
        seq<-FASTA[start_s:end_s]
      }else if(start_s<0 && end_s<length(FASTA)){
        start_s<-length(FASTA)+start_s
        seq<-c(FASTA[start_s:length(FASTA)],FASTA[1:end_s])
      }else if(start_s>0 && end_s>length(FASTA)){
        end_s<-end_s-length(FASTA)
        seq<-c(FASTA[start_s:length(FASTA)],FASTA[1:end_s])
      }
      seq<-paste(seq,collapse="")
    
      Sequence[i]<-seq
    }
    finaltab<-cbind(finaltab,Sequence)
  }
  
  
  return(finaltab)
}

# Function to merge microarray and RNAseq data, normalize with Combat and create final table 
interGSE_analysis<-function(path_input,tmp_MA_path,path_rna,results_path,path_fulltab_res,path_output,path_operons,tmp_genome_path,OPERONS_DETECTION,asm_pref,multispecies){
  
  # Merge microarray and RNAseq data
output<-unique_matrix(tmp_MA_path, path_rna, tmp_genome_path)
  if(length(output)>0){
    tabella_completa<-output[[1]]
    batches<-output[[2]]
    unique_batches<-unique(batches)
    
    ind_col<-which(regexpr("GSE",colnames(tabella_completa))>0)
    E<-tabella_completa[,ind_col]
    ESP<-matrix(as.numeric(E),nrow=nrow(E),ncol=ncol(E))
    colnames(ESP)<-colnames(E)
    rm(E)
    
    # Divide experiment according to specie
    file_input<-dir("/input")
    GSE_tab_complete<-c()
    for(i in 1:length(file_input)){
      text<-read.xlsx(file.path("/input",file_input[i]))
      GSE_tab<-cbind(text$Accession,text$taxon)
      if(length(GSE_tab_complete)==0){
        GSE_tab_complete<-GSE_tab
      }else{
        GSE_tab_complete<-rbind(GSE_tab_complete,GSE_tab)
      }
    }

    if(multispecies==TRUE){
       unique_species<-unique(GSE_tab_complete[,2])
    
      for(i in 1:length(unique_species)){
        ind<-which(unique_species[i]==GSE_tab_complete[,2])

        for(j in 1:length(ind)){
          ind1<-which(regexpr(GSE_tab_complete[ind[j],1],colnames(ESP))>0)
          if(length(ind1)>0){
            if(j==1){
              ESP1<-ESP[,ind1]
              batches1<-rep(GSE_tab_complete[ind[j],1],length(ind1))
            }else{
              ESP1<-cbind(ESP1,ESP[,ind1])
              batches1<-c(batches1,rep(GSE_tab_complete[ind[j],1],length(ind1)))
            }
          }
        }


        # INTRAEXPERIMENTAL NORMALIZATION: COMBAT
        if(length(ind)>1){
          combat_edata1 = ComBat(dat=ESP1, batch=batches1, mod=NULL, par.prior=T, prior.plots=F, mean.only = T)
        }else if(length(ind)==1){
          combat_edata1=ESP1
        }
        
        if(i==1){
          final_ESP<-combat_edata1
        }else{
          final_ESP<-cbind(final_ESP,combat_edata1)
        }
      }
    }else{
      batches1<-c()
      for(i in 1:nrow(GSE_tab_complete)){
        ind1<-which(regexpr(GSE_tab_complete[i,1],colnames(ESP))>0)
        batches1<-c(batches1,rep(GSE_tab_complete[i,1],length(ind1)))
      }
      combat_edata1 = ComBat(dat=ESP, batch=batches1, mod=NULL, par.prior=T, prior.plots=F, mean.only = T)
      final_ESP<-combat_edata1
    }

    # Final table 
    finaltab<-create_final_table(final_ESP,tabella_completa)
    
    # ADD promoter sequence
    finaltab<-add_promoter_sequence(path_output,finaltab,asm_pref)
    
    write.xlsx(finaltab,file.path(results_path,"final_table.xlsx"))
    
    # copy of full tabs file and remove of tmp directory
    last_operations(tmp_MA_path,path_rna,path_fulltab_res,path_output)
  }else{
    message("######## The analysis ended with an error. ########")
  }
}
  
