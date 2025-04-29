# Reprocessing of OrthoFinder output
find_ortologhi<-function(geni_ortologhi,pref,other,asm_pref,asm_code,tabella_finale,n_amin,threashold){
  
  # pref
  ind_pref<-which(asm_pref==colnames(geni_ortologhi))
  
  # other
  ind_other<-which(gsub("-","\\.",asm_code)==colnames(geni_ortologhi))
  if(is_empty(ind_other)){
    ind_other<-which(gsub(" ","_",asm_code)==colnames(geni_ortologhi))
  }
  
  # remove empty orthogroups
  ind<-c(which(geni_ortologhi[,ind_pref]==""),which(geni_ortologhi[,ind_other]==""))
  if(length(ind)>0){
    geni_ortologhi<-geni_ortologhi[-ind,]
  }
  
  # Considering only not empty orthogroups
  if(length(geni_ortologhi)>0 && nrow(geni_ortologhi)>0){
    # Find orthogroups containing more than two genes
    ind1<-which(regexpr(",",geni_ortologhi[,ind_pref])>0)
    ind2<-which(regexpr(",",geni_ortologhi[,ind_other])>0)
    ind<-c(ind1,ind2)
    if(length(ind)>0){
      tab<-geni_ortologhi[ind,] # > 2 genes
      tab1_1<-geni_ortologhi[-ind,] # = 2 genes
      
      gene_v<-c()
      elim<-c()
      for(i in 1:nrow(tab)){
        g_1<-strsplit(tab[i,ind_pref],", ")[[1]]
        l_1<-length(g_1)
        
        g_2<-strsplit(tab[i,ind_other],", ")[[1]]
        l_2<-length(g_2)
        
        # The alignment scores for each protein belonging to an orthogroup are 
        # obtained, and the proteins with the highest scores and a difference in 
        # the number of amino acids < n_amin are associated.
        score<-matrix(nrow=l_1,ncol=l_2)
        score_v<-c()
        similarity_v<-c()
        similarity<-matrix(nrow=l_1,ncol=l_2)
        for(j in 1:l_1){
          i_1<-which(g_1[j]==pref[,2])
          for(z in 1:length(i_1)){
            prot1<-pref[i_1[z],5]
            for(k in 1:l_2){
              i_2<-which(g_2[k]==other[,2])
              for(x in 1:length(i_2)){
                prot2<-other[i_2[x],5]
                if(abs(nchar(prot1)-nchar(prot2))<=n_amin){
                  globalAlign <-pairwiseAlignment(prot1,prot2, type = "global")
                  matches <- nmatch(globalAlign)               
                  aligned_length <- nchar(pattern(globalAlign))  
                  if(length(score_v)==0){
                    score_v<-globalAlign@score
                    similarity_v <- (matches / aligned_length) * 100
                  }else{
                    score_v<-c(score_v,globalAlign@score)
                    similarity_v<-c(similarity_v,(matches / aligned_length) * 100)
                  }
                }else{
                  if(length(score_v)==0){
                    score_v<-0
                    similarity_v<-0
                  }else{
                    score_v<-c(score_v,0)
                    similarity_v<-c(similarity_v,0)
                  }
                }
              }
              score[j,k]<-max(score_v)
              ind_sim<-which(score_v==max(score_v))
              similarity[j,k]<-similarity_v[ind_sim]
              score_v<-c()
              similarity_v<-c()
            }
          }
        }
        
        while(length(score)>1 || l_1>1 || l_2>1){
          if(max(score)==0){
            elim<-c(elim,i)
            break
          }
          ind<-which(score==max(score),arr.ind = T)
          if(similarity[ind[1,1],ind[1,2]]>=threashold){
            gene1<-g_1[ind[1,1]]
            gene2<-g_2[ind[1,2]]
            g_1<-g_1[-ind[1,1]]
            g_2<-g_2[-ind[1,2]]
            l_1<-l_1-1
            l_2<-l_2-1
            if(ind_pref>ind_other){
              gene_o<-cbind(gene2,gene1)
            }else{
              gene_o<-cbind(gene1,gene2)
            }
            # final association
            gene_v<-rbind(gene_v,gene_o)
            if(l_1>1 && l_2>1){
              score<-score[-ind[1],]
              similarity<-similarity[-ind[1],]
            }else if(l_1>1 && l_2==1){
              score<-score[-ind[1]]
              similarity<-similarity[-ind[1]]
            }else{
              break
            }
            if(l_2>1 && l_1>1){
              score<-score[,-ind[2]]
              similarity<-similarity[,-ind[2]]
            }else if(l_2>1 && l_1==1){
              score<-score[-ind[2]]
              similarity<-similarity[-ind[2]]
            }else{
              break
            }
          }else{
            elim<-c(elim,i)
            break
          }
        }
      }
      if(length(elim)>0){
        tab<-tab[-elim,]
      }
      ind<-c()
      for(i in 1:nrow(gene_v)){
        ind[i]<-which(regexpr(gene_v[i,1],tab[,4])>0)
      }
      tab<-cbind(tab[ind,1:3],gene_v)
      colnames(tab)<-colnames(tab1_1)
    }else{
      tab1_1<-geni_ortologhi
      tab<-c()
    }
    
    # orthogroups containing only two genes
    if(nrow(tab1_1)>0){
      index<-c()
      dist<-c()
      for(i in 1:nrow(tab1_1)){
        i_1<-which(tab1_1[i,ind_pref]==pref[,2])
        i_2<-which(tab1_1[i,ind_other]==other[,2])
        dist<-c()
        similarity_v<-c()
        for(j in 1:length(i_1)){
          prot1<-pref[i_1[j],5]
          for(k in 1:length(i_2)){
            prot2<-other[i_2[k],5]
            if(length(dist)==0){
              dist<-abs(nchar(prot1)-nchar(prot2))
            }else{
              dist<-c(dist,abs(nchar(prot1)-nchar(prot2)))
            }
            globalAlign <-pairwiseAlignment(prot1,prot2, type = "global")
            matches <- nmatch(globalAlign)               
            aligned_length <- nchar(pattern(globalAlign))  
            if(length(similarity_v)==0){
              similarity_v <-(matches / aligned_length) * 100
            }else{
              similarity_v<-c(similarity_v,(matches / aligned_length) * 100)
            }
          }
        }
        dist<-min(dist)
        similarity_v<-max(similarity_v)
        if(dist>=n_amin || similarity_v<threashold){
          index<-c(index,i)
        }
      }
      if(length(index)>0){
        tab1_1<-tab1_1[-index,]
      }
    }
    if (nrow(tab1_1)>0 && nrow(tab)){
      # Final couples of orthologous genes
      geni_ortologhi<-rbind(tab1_1,tab)
    }else if(nrow(tab1_1)>0 && nrow(tab)==0){
      geni_ortologhi<-tab1_1
    }else if(nrow(tab1_1)==0 && nrow(tab)>0){
      geni_ortologhi<-tab
    }
    
    if(nrow(geni_ortologhi)>0){
      # Create final table
      new<-matrix(ncol=3,nrow = nrow(tabella_finale))
      for(i in 1:nrow(geni_ortologhi)){
        ind<-which(geni_ortologhi[i,ind_other]==other[,2])
        riga_other<-other[ind,2:4]
        ind<-which(geni_ortologhi[i,ind_pref]==tabella_finale[,2])
        new[ind,1:length(as.character(riga_other))]<-as.character(riga_other)
      }
      tabella_finale<-cbind(tabella_finale,new)
    }
    
  }
  
    
  return(tabella_finale)
}

# Create a file in .fa format associating gene and proteins
create_genome_table<-function(asm_gb,path_file,orthopath,taxid_species){
  # read gb file
  SS<-readLines(asm_gb)
  genome<-read_gbfile(SS)
  
  file<-c()
  for(i in 1:nrow(genome)){
    if(i==1 && is.na(genome[i,5])==0){
      file<-rbind(paste(">",genome[i,2],sep=""),genome[i,5])
    }else if(is.na(genome[i,5])==0){
      file<-rbind(file,paste(">",genome[i,2],sep=""),genome[i,5])
    }
  }
  asm<-gsub(".gb","",asm_gb)
  asm<-gsub(" ","_",asm)
  name_file<-paste(asm,".fa",sep="")
  file<-data.frame(file)
  fwrite(file,name_file,col.names = F)
  
  # Convert path in Ubuntu format
  setwd(path_file)
  path_file1<-gsub("[A-Za-z ]+:/[A-Za-z ]+/[A-Za-z ]+","..",path_file)
  path2<-gsub("[A-Za-z_0-9]+","..",path_file1)
  
  # Create a folder to copy the .fa files that Orthofinder will process.
  testo0<-paste('bash -c',shQuote(paste("mkdir -p ","../../../",path2, orthopath, "/ExampleData/", taxid_species,"/",sep="")),sep=" ")
  system(testo0)
  testo0<-paste('bash -c',shQuote(paste("cp ", name_file," ../../../",path2, orthopath, "/ExampleData/", taxid_species,"/",sep="")),sep=" ")
  system(testo0)

  # Remove file from the current directory
  file.remove(paste(path_file,name_file,sep="/"))

  return(genome)
}

# Function that extracts important information from the .gb file
read_gbfile<-function(SS){
  
  strings<-c("CDS + [0-9complement]+."," +gene + [<0-9complement]+.+","/locus_tag=","/old_locus_tag=","/translation=","/gene=","/gene_synonym=","join","ORIGIN")
  
  # Identify the lines containing the strings "CDS" and "gene," and remove the 
  # lines related to "join." Restrict the search for gene information to the 
  # lines between a "CDS" and the next "gene."
  ind_cds<-which(regexpr(strings[1],SS)>0)
  ind_gene<-which(regexpr(strings[2],SS)>0)
  ind_join<-which(regexpr(strings[8],SS[ind_cds])>0)
  if(length(ind_join)>0){
    ind_cds<-ind_cds[-ind_join]
  }
  ind_join<-which(regexpr(strings[8],SS[ind_gene])>0)
  if(length(ind_join)>0){
    ind_gene<-ind_gene[-ind_join]
  }
  ind_gene<-c(ind_gene,which(regexpr(strings[9],SS)>0))
  
  limits<-matrix(nrow = length(ind_cds),ncol = 2)
  for(i in 1:nrow(limits)){
    limits[i,1]<-ind_cds[i]
    limits[i,2]<-ind_gene[min(which(ind_gene>ind_cds[i]))]
  }

  gene_names<-c()
  locus_tag<-c()
  old_locus_tag<-c()
  synonim<-c()
  proteins<-c()
  for(i in 1:nrow(limits)){
    lines<-c()
    lines<-SS[limits[i,1]:limits[i,2]]
    # gene name
    ind1<-which(regexpr(strings[6],lines)>0)[1]
    if(length(ind1)==0){
      gene_names[i]<-NA
    }else{
      a<-gsub(strings[6],"",lines[ind1])
      gene_names[i]<-gsub("[^a-zA-Z0-9]","",a)
    }
    # locus tag
    ind2<-which(regexpr(strings[3],lines)>0)[1]
    if(length(ind2)==0){
      locus_tag[i]<-NA
    }else{
      b<-gsub(strings[3],"",lines[ind2])
      locus_tag[i]<-gsub("[^a-zA-Z0-9_]","",b)
    }
    # old locus tag
    ind3<-which(regexpr(strings[4],lines)>0)[1]
    if(length(ind3)==0){
      old_locus_tag[i]<-NA
    }else{
      c<-gsub(strings[4],"",lines[ind3])
      old_locus_tag[i]<-gsub("[^a-zA-Z0-9_]","",c)
    }
    # synonims
    ind4<-which(regexpr(strings[7],lines)>0)[1]
    if(length(ind4)==0){
      synonim[i]<-NA
    }else{
      d<-gsub(strings[7],"",lines[ind4])
      synonim[i]<-gsub("[^a-zA-Z0-9;]","",d)
    }
    # proteina
    ind5<-which(regexpr(strings[5],lines)>0)[1]
    if(length(ind5)==0 || is.na(ind5)){
      proteins[i]<-NA
    }else{
      ind6<-which(regexpr("\"",lines)>0)
      ind6<-ind6[which(ind6>=ind5)]
      if(length(ind6)>1 || length(ind6)==0){
        ind6<-min(ind6)
      }
      e<-gsub(strings[5],"",lines[ind5:ind6])
      e<-paste0(e,collapse = "")
      proteins[i]<-gsub("[^a-zA-Z0-9]","",e)
    }
  }
  gen<-cbind(gene_names,locus_tag,old_locus_tag,synonim,proteins)
  colnames(gen)<-c("gene","locus_tag","old_locus_tag","synonim","sequence_AA")
  return(gen)
}

# Function to find orthologous genes
main_ortologhi<-function(path_input,path_output,asm_pref,tmp_genome_path,orthopath, THREADS,n_amin,threashold){
  source("/scripts/Microarray/dowload_file_assembly.R")
  
  setwd(path_input)
  file_input<-dir(path_input)
  asm_code<-c()
  for(i in 1:length(file_input)){
    text<-read.xlsx(file_input[i])
    asm_code<-c(asm_code,text$reference_genome)
  }
  asm_code<-unique(asm_code)
  
  message("..downloading .gb files for each strain..")
  # download of gb and fasta files for each strain by ftp address
  reading_testo(path_output,tmp_genome_path,asm_code,asm_pref)
  
  
  # Directory to save the input of OrthoFinder 
  name_dir<-"ortho_output"
  # gb file
  asm_gb<-dir(tmp_genome_path,pattern=".gb")
  
  setwd(tmp_genome_path)
  # preferred genome
  asm_pref_gb<-asm_gb[which(regexpr(asm_pref,asm_gb)>0)]
  pref<-create_genome_table(asm_pref_gb,tmp_genome_path,orthopath,name_dir)
  ind<-which(asm_pref_gb==asm_gb)
  asm_gb<-asm_gb[-ind]
  # asm code of other genomes
  asm_code<-gsub(".gb","",asm_gb)
  
  # Initialize the final table with the information of the preferred genome
  locus<-unique(pref[,2])
  tabella_finale<-matrix(ncol = 4,nrow = length(locus))
  for(i in 1:length(locus)){
    ind<-which(pref[,2]==locus[i])[1]
    tabella_finale[i,]<-pref[ind,1:4]
  }

  for(i in 1:length(asm_gb)){

    other<-create_genome_table(asm_gb[i],tmp_genome_path,orthopath,name_dir)
    
    # Directory of results
    d<-date()
    month<-gsub(" ","",substr(d,5,7))
    day<-gsub(" ","",substr(d,9,10))
    if(nchar(day)==1){
      day<-paste("0",day,sep="")
    }
    results_dir<-paste("Results_",month,day,sep = "")
    
    message(paste("Analyzed strain:",asm_code[i],sep=" "))
    
    # Create the directory of the input of OrthoFinder
    testo0<-paste('bash -c',shQuote(paste("mkdir -p ",orthopath,"/ExampleData/",name_dir,'/',sep="")),sep=" ")
    system(testo0)
    # Run OrthoFinder
    testo0<-paste('bash -c',shQuote(paste("cd ",orthopath," && ./orthofinder -t ", THREADS, " -f ExampleData/",name_dir,"/",sep="")),sep=" ")
    system(testo0)
    # Copy the results of orthoFinder in the current directory
    path_output<-paste(tmp_genome_path,"output",sep="/")
    if(dir.exists(path_output)==0){
      dir.create(path_output)
    }
    testo0 <- paste('bash -c', shQuote(paste("cd", paste0(orthopath, "/ExampleData/", name_dir, "/OrthoFinder/", results_dir, "/Phylogenetic_Hierarchical_Orthogroups"), "&& cp N0.tsv", path_output)))
    system(testo0)
    # Remove the directory containing the output of OrthoFinder
    testo0<-paste('bash -c',shQuote(paste("cd ",orthopath,"/ExampleData/",name_dir,'/OrthoFinder/ ',"&& rm -r ",results_dir,sep="")),sep=" ")
    system(testo0)
    # Remove the .fa file of the "other" genome from the input directory
    file<-paste(gsub(" ","_",asm_code[i]),".fa",sep="")
    testo0<-paste('bash -c',shQuote(paste("cd ",orthopath,"/ExampleData/",name_dir,"&& rm ",file,sep="")),sep=" ")
    system(testo0)
    
    Orthofile<-paste(path_output,"/N0.tsv",sep = "")
    geni_ortologhi<-read.delim(Orthofile,sep = '\t', header = TRUE, fill = T)
    # Final table
    
    tabella_finale<-find_ortologhi(geni_ortologhi,pref,other,asm_pref,asm_code[i],tabella_finale,n_amin,threashold)

  }
  tabella_finale<-remove_empty(tabella_finale,"cols")
  # Write the final reprocessed output 
  name_file<-"orthologous_genes.txt"
  write.table(tabella_finale,paste("/output",name_file,sep="/"),sep="\t",quote = F,col.names = F,row.names = F)
  
  #return(geni_ortologhi_tot)
}
