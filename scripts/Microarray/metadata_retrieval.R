# Function that retrieves data related to the analyzed experiments
data_retrieval<-function(GSE,path_file,asm_name,taxid_species){
  
  n_GSE<-length(GSE)
  
  GPL<-c()
  GSM<-matrix(ncol = n_GSE,nrow = 1000)
  n_GSM<-c()
  n_canali<-c()
  manufacturer<-c()
  organism<-c()
  taxid<-c()
  
  locus_tag<-matrix(ncol = n_GSE,nrow = 50000)
  gene_name<-matrix(ncol = n_GSE,nrow = 50000)
  nomi_sonda<-matrix(ncol = n_GSE,nrow = 50000)
  
  for(j in 1:n_GSE){
    a<-getGEO(GSE[j],AnnotGPL = T,GSEMatrix = F)
    # GPL
    GPL[j]<-a@header$platform_id
    # GSM
    GSM_code<-c()
    GSM_code<-a@header$sample_id
    n_GSM[j]<-length(GSM_code)
    # GSM in each GSE
    GSM[1:length(GSM_code),j]<-GSM_code
    # channels
    n_canali[j]<-a@gsms[[1]]@header$channel_count
    # manufacturer
    manufacturer[j]<-a@gpls[[1]]@header$manufacturer
    # organism
    organism[j]<-a@gpls[[1]]@header$organism
    # taxid
    taxid[j]<-a@gpls[[1]]@header$taxid
    # attributes GPL
    X<-colnames(a@gpls[[1]]@dataTable@table)
    if(length(X)>0){
      # locus tag
      for (i in 1:length(X)){
        if(X[i]=="ORF_1"){
          locus_tag[1:length(a@gpls[[1]]@dataTable@table[,i]),j]<-a@gpls[[1]]@dataTable@table[,i]
          break
        }else if(i!=length(X) & X[i+1]!="ORF_1" & (X[i]=="ORF" || X[i]=="NAME")){
          locus_tag[1:length(a@gpls[[1]]@dataTable@table[,i]),j]<-a@gpls[[1]]@dataTable@table[,i]
          break
        }else if(X[i]=="ORF" || X[i]=="NAME" || X[i]=="gene_name" || X[i]=="SystematicName" || X[i]=="ORF_LIST"){
          locus_tag[1:length(a@gpls[[1]]@dataTable@table[,i]),j]<-a@gpls[[1]]@dataTable@table[,i]
          break
        } else if(regexpr("LOCUS_[A-Z]+[0-9]+",X[i])==1){
          locus_tag[1:length(a@gpls[[1]]@dataTable@table[,i]),j]<-a@gpls[[1]]@dataTable@table[,i]
          break
        }
      }
      # gene name
      if(length(which(X=="Gene Symbol"))){
        gene_name[1:length(a@gpls[[1]]@dataTable@table$`Gene Symbol`),j]<-a@gpls[[1]]@dataTable@table$`Gene Symbol`
      }else{
        gene_name[1:length(locus_tag[,j]),j]<-locus_tag[,j]
      }
      # probe name
      if (length(which(X=="ID"))>0){
        nomi_sonda[1:length(a@gpls[[1]]@dataTable@table$ID),j]<-a@gpls[[1]]@dataTable@table$ID
      }else{
        nomi_sonda[1:length(locus_tag[,j]),j]<-locus_tag[,j]
      }
    }
  }
  nomi_sonda<-remove_empty(nomi_sonda,"rows")
  locus_tag<-remove_empty(locus_tag,"rows")
  gene_name<-remove_empty(gene_name,"rows")
  
  # Remove GSE missing information
  ind<-c()
  for(i in 1:ncol(locus_tag)){
    if(length(unique(locus_tag[,i]))==1 && is.na(unique(locus_tag[,i]))==1){
      ind[i]<-1
    }
  }
  ind<-which(ind==1)
  if(length(ind)>0){
    locus_tag<-locus_tag[,-ind]
    gene_name<-gene_name[,-ind]
    GPL<-GPL[-ind]
    GSE<-GSE[-ind]
    manufacturer<-manufacturer[-ind]
    n_canali<-n_canali[-ind]
    organism<-organism[-ind]
    taxid<-taxid[-ind]
    n_GSE<-length(GSE)
    nomi_sonda<-nomi_sonda[,-ind]
    GSM<-GSM[,-ind]
    n_GSM<-n_GSM[-ind]
    asm_name<-asm_name[-ind]
    taxid_species<-taxid_species[-ind]
  }
  
  # Unique GPL
  GPL_uniche<-unique(GPL)
  n_GPL<-length(GPL_uniche)
  # For each GPL, create a file reporting ID_probe, locus tag and gene name
  for(i in 1:n_GPL){
    ind<-which(GPL==GPL_uniche[i])
    if(length(ind)>1){
      ind<-ind[1]
    }
    tab<-cbind(nomi_sonda[,ind],locus_tag[,ind],gsub("[^A-Za-z1-9_/ ]","",gene_name[,ind]))
    tab<-remove_empty(tab,"rows")
    colnames(tab)<-c("ID_probe","locus","gene_name")
    tab<-data.frame(tab)
    write.table(tab,paste(path_file,"/",GPL_uniche[i],".txt",sep=""),sep="\t",quote=F,row.names = F)
  }
  
  tabella1<-matrix(ncol=n_GSE,nrow = 100000)
  # Download data for each GSE
  for(i in 1:n_GSE){
    list_file<-c()
    list_file_new<-c()
    list_file_new2<-c()
    path1<-c()
    path2<-c()
    directory<-c()
    ind<-c()
    
    if(dir.exists(paste(path_file,GSE[i],sep="/"))==0 || length(dir(paste(path_file,GSE[i],sep="/"),pattern = "GSM"))!=n_GSM[i]){
      
      for(ii in 1:3){
        tryCatch({
          b<-getGEOSuppFiles(GSE[i],makeDirectory = TRUE,baseDir = path_file)
          path1 <- paste(path_file, GSE[i], sep = "/")
          directory <- dir(path1)
          ind <- which(regexpr("_RAW.tar",directory)<0)
          if(length(ind) > 0) {
            file.remove(paste(path1,directory[ind], sep = "/"))
            directory <- directory[-ind]
            unlink(paste(path1, "immagini", sep = "/"),force = TRUE, recursive = T)
          }
          break
        }, error = function(e) {
          message("An error occurred while downloading file for",GSE[i])
          Sys.sleep(60)
          
        })
      }
      if (ii == 3) {
        next
      }
      if(length(directory)>0){
        path2<-paste(path1,directory,sep="/")
        untar(path2,exdir=path1)
        
        exprreg<-"^GPL[0-9]+"
        list_file<-dir(path1)
        ind<-which(regexpr(exprreg,list_file)>0)
        if(length(ind)>0){
          list_file_new2<-gsub(".gz","",list_file[ind])
          for(ss in 1:length(ind)){
            gunzip(filename=paste(path1,list_file[ss],sep="/"),destname = paste(path1,list_file_new2[ss],sep="/"),overwrite=T)
            tab<-fread(paste(path1,list_file_new2[ss],sep="/"))
            nomi_colonne<-colnames(tab)
            if(length(which("GeneSymbol"==nomi_colonne))>0){
              tabella1[1:nrow(tab),i]<-tab$GeneSymbol
            }
          } 
          list_file<-list_file[-ind]
        }
        
        if(regexpr("[Aa]ffymetrix",manufacturer[i])>0){
          list_file<-dir(path1)
          expr<-"^GSM.+.CEL.+$"
          ind<-which(regexpr(expr,list_file)<0)
        }else{
          list_file<-dir(path1)
          expr<-"^GSM.+.pair.+$"
          ind<-which(regexpr(expr,list_file)<0)
          if(length(ind)==length(list_file)){
            expr<-"^GSM.+.txt.+$"
            ind<-which(regexpr(expr,list_file)<0)
          }
          if(length(ind)==length(list_file)){
            expr<-"^GSM.+.csv.+$"
            ind<-which(regexpr(expr,list_file)<0)
          }
          if(length(ind)==length(list_file)){
            expr<-"^GSM.+.gpr.+$"
            ind<-which(regexpr(expr,list_file)<0)
          }
        }
        expr<-"RMA"
        ind1<-which(regexpr(expr,list_file)>0)
        ind<-c(ind,ind1)
        file.remove(paste(path1,list_file[ind],sep="/"))
        list_file_new<-dir(path1,".gz")
        list_file_new2<-gsub(".gz","",list_file_new)
        for(k in 1:length(list_file_new)){
          gunzip(filename=paste(path1,list_file_new[k],sep="/"),destname = paste(path1,list_file_new2[k],sep="/"),overwrite=T)
        }
      }else{
        if(dir.exists(path1)){
          testo0<-paste('bash -c',shQuote(paste("rmdir ",path1, sep="")),sep=" ")
          system(testo0)
        }else{
          message(paste("#### No supplementary files are available for ",GSE[i]," ####"))
        }
        
      }
    }
  }
  
  # Get a list of files
  files <- list.files(path_file, pattern = "^GSE")
  GSE_rimasti <- grep("full_tab", files, value = TRUE, invert = TRUE)
  ind <- numeric(length(GSE_rimasti))
  for (i in seq_along(GSE_rimasti)) {
    current_index <- which(GSE == GSE_rimasti[i])
    if (length(current_index) > 0) {
      ind[i] <- current_index[1]
    }
  }

  if(length(ind)>0){
    GPL<-GPL[ind]
    GPL_uniche<-unique(GPL)
    if(is.matrix(GSM)){
      GSM<-GSM[,ind]
    }    
    n_canali<-n_canali[ind]
    manufacturer<-manufacturer[ind]
    organism<-organism[ind]
    taxid<-taxid[ind]
    if(is.matrix(tabella1)){
      tabella1<-tabella1[,ind]
    }
    
    GSE<-GSE_rimasti
    n_GSM<-n_GSM[ind]
    if(is.matrix(locus_tag)){
      locus_tag<-locus_tag[,ind]
    }
    if(is.matrix(gene_name)){
      gene_name<-gene_name[,ind]
    }
    if(is.matrix(nomi_sonda)){
      nomi_sonda<-nomi_sonda[,ind]
    }
    
    asm_name<-asm_name[ind]
    taxid_species<-taxid_species[ind]

    # metadata information
    tab1<-data.frame(GSE,GPL,organism,taxid,taxid_species,asm_name,n_canali,n_GSM,manufacturer)
    write.xlsx(tab1,paste(path_file,"/total_metadata_GSE.xlsx",sep=""))
  
    # locus_tag
    if(is.matrix(locus_tag)){
      locus_tag<-remove_empty(locus_tag,"rows")
    }else{
      ind<-which(is.na(locus_tag))
      if(length(ind)>0){
        locus_tag<-locus_tag[-ind]
      }
      
    }
    tab2<-data.frame(locus_tag)
    colnames(tab2)<-GSE
    write.xlsx(tab2,paste(path_file,"/total_locus_tag.xlsx",sep=""))
    
    # gene name
    if(is.matrix(gene_name)){
      gene_name<-remove_empty(gene_name,"rows")
    }else{
      ind<-which(is.na(gene_name))
      if(length(ind)>0){
        gene_name<-gene_name[-ind]
      }
    }
    tab3<-data.frame(gene_name)
    colnames(tab3)<-GSE
    write.xlsx(tab3,paste(path_file,"/total_gene_name.xlsx",sep=""))
    
    # GSM
    if(is.matrix(GSM)){
      GSM<-remove_empty(GSM,"rows")
    }else{
      ind<-which(is.na(GSM))
      if(length(ind)>0){
        GSM<-GSM[-ind]
      }
      
    }
    tab4<-data.frame(GSM)
    colnames(tab4)<-GSE
    write.xlsx(tab4,paste(path_file,"/total_GSM.xlsx",sep=""))
    
    # probe name
    if(is.matrix(tabella1)){
      tabella1<-remove_empty(tabella1,"rows")
    }else{
      ind<-which(is.na(tabella1))
      if(length(ind)>0){
        tabella1<-tabella1[-ind]
      }
      
    }
    tab5<-data.frame(tabella1)
    colnames(tab5)<-GSE
    write.xlsx(tab5,paste(path_file,"/total_probes.xlsx",sep=""))

  }
  
  
}

# Retrieving experimental conditions for each GSE
Groups_retrieval = function(gsm,GSE,path_GSE){ 
  
  cond= list()
  # Retrieving experimental conditions for each GSM
  for (i in 1:length(gsm)){
    for(j in 1:3){
      tryCatch({
      gsmsingolo=getGEO(GEO= gsm[i], destdir=path_GSE)
        if(length(gsmsingolo)>0){
          break
        }
      },error=function(e){
          message("An error occurred while downloading file for",gsm[i])
          Sys.sleep(60)
          
      })
    }
    if(j==3){
      next
    }
    cond[[i]]= paste(gsmsingolo@header$title, gsmsingolo@header$characteristics_ch1)
    cond[[i]]<-gsub("[;:/%*()]"," ",cond[[i]])
    
    
    cond[[i]]= str_replace_all(cond[[i]], "\\+", " plus ")
    cond[[i]]= str_replace_all(cond[[i]], "\\-", " minus ")
    cond[[i]]= paste(cond[[i]], collapse='')
  }
  soft_file<-dir(path_GSE,pattern=".soft")
  file.remove(paste0(path_GSE, "/",soft_file))
  cond<- unlist(cond)
  
  # Computing distances between strings to cluster experimental conditions within the same group
  if(length(cond)>1){
    distances= adist(cond)
    rownames(distances) <- colnames(distances) <- cond
    threshold<- min(distances[distances > 0])+1 
    
    hc <- hclust(as.dist(distances))
    groups<- data.frame(cond,cutree(hc,h=threshold)) 
    
  } else { groups= data.frame(cond[[1]], '1') }
  
  colnames(groups)<- c("Condition", "Group")
  lis<-list(groups)
  return(lis)
}

# Function that generates a file with the experimental conditions for each GSE and GSM
do_file_for_group_control<-function(path_file){
  
  tab1<-read.xlsx(file.path(path_file,"total_metadata_GSE.xlsx"))
  GSE<-tab1$GSE
  GPL<-tab1$GPL
  organism<-tab1$organism
  taxid<-tab1$taxid
  n_canali<-tab1$n_canali
  n_GSM<-tab1$n_GSM
  manufacturer<-tab1$manufacturer
  asm_name<-tab1$asm_name
  
  GSM<-read.xlsx(file.path(path_file,"total_GSM.xlsx"))
  
  tabella<-matrix(ncol =7 ,nrow = length(GSE))
  
  n_GSM_tot<-c()
  insieme<-c()
  for(i in 1:ncol(GSM)){
    path_GSE<-file.path(path_file,GSE[i])
    #setwd(path_GSE)
    gsm<-GSM[1:n_GSM[i],i]
    
    GSM_cond<-Groups_retrieval(gsm,GSE[i],path_GSE)[[1]]
    accession.number<-paste(gsm,collapse=",")
    groups<-paste(GSM_cond$Group,collapse = ",")
    Experimental_condition<-paste(GSM_cond$Condition,collapse = ",")
    GSEcode<-GSE[i]
    asm_namev<-asm_name[i]
    taxidv<-taxid[i]
    m.o.name<-organism[i]
    tabella[i,]<-cbind(m.o.name,taxidv,asm_namev,GSEcode,Experimental_condition,groups,accession.number)
    n_GSM_tot[i]<-length(gsm)
  }
  tabella<-data.frame(tabella)
  colnames(tabella)<-c("m.o.name",	"taxid",	"asm_name",	"GSE",	"Experimental_condition",	"groups",	"accession.number")
  
  file_name<-paste(path_file,"gsm_condition.xlsx",sep="/")
  write.xlsx(tabella,file_name)

}

# Function that downloads all files related to gene expression in the analyzed experiments
metadata_retrieval<-function(path_file,file_MA,path_input,path_output){
 # setwd(path_file)
  
  GSE<-c()
  asm_name<-c()
  if(length(file_MA)>1){
    for(i in 1:length(file_MA)){
      text<-read.xlsx(paste("/input/",file_MA[i],sep=""))
      GSE<-c(GSE,text$Accession)
      asm_name<-c(asm_name,text$reference_genome)
    }
  }else{
    text<-read.xlsx(paste("/input/",file_MA,sep=""))
    GSE<-c(text$Accession)
    asm_name<-c(text$reference_genome)
  }
  
  
  assembly_summary<-fread(file.path(path_output,"assembly_summary.txt"), quote=F)
  taxid_species<-c()
  for(i in 1:length(asm_name)){
    ind<-which(assembly_summary$asm_name==asm_name[i])
    taxid_species[i]<-assembly_summary$species_taxid[ind]
  }
  n_GSE<-length(GSE)
  
  # Retrieve of metadata and expression files
  risultati1<-data_retrieval(GSE,path_file,asm_name,taxid_species)
  # Create group file 
  if(!file.exists("/output/output/tmp/microarray/gsm_condition.xlsx")){
    risultati2<-do_file_for_group_control(path_file)
  }else{
    message("Success: Experimental Groups already defined!")
  }
  
  
  
  cat('Before proceeding, check if the groups and replicas have been assigned correctly within file "gsm_condition.txt". 
If not, please edit the file and overwrite it. When you are done, press enter:')
  line <- readline()
}

