# Function to download file fna -->  .gbff
download<-function(ftp_path,path_dir_org,asm){
  
    ftp_string<- gsub(".*/", "", ftp_path)
    GFF_ftp<- paste0(ftp_path, "/", ftp_string, "_genomic.gbff.gz")
    
    annot<- curl_download(GFF_ftp, destfile=paste0(path_dir_org,"/", ftp_string, "_genomic.gbff.gz"), mode="wb")
    new_name<-gsub(".+gbff.gz",paste(path_dir_org,"/",asm,".gb",sep=""),annot)
    gunzip(filename=annot,destname=new_name,overwrite=T)
 
}

# Function to read assembly_summary and to download file for each version of genome  
reading_testo<-function(path_output,tmp_genome_path,asm_code,asm_pref){

  if (!file.exists(file.path("/output/assembly_summary.txt"))){
      message("Downloading assembly summary from RefSeq..")
      #summary<- curl_download("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt" , destfile=file.path("/output/assembly_summary.txt"), mode="wb")
      summary<- curl_download("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt" , destfile=file.path(path_output,"/assembly_summary.txt"), mode="wb")
  }

  #summary<- fread(file.path("/output/assembly_summary.txt"), quote=F)
  summary<- fread(file.path(path_output,"assembly_summary.txt"), quote=F)
  
  ind<-c()
  taxid_species<-c()
  for(i in 1:length(asm_code)){
    ind<-which(summary$asm_name==asm_code[i])
    taxid_species[i]<-summary$species_taxid[ind]
  }
  ind_pref<-which(summary$asm_name==asm_pref)
  taxid_species<-unique(summary$species_taxid[ind],summary$species_taxid[ind_pref])
  
  ind<-c()
  for(i in 1:length(taxid_species)){
    ind<-c(ind,which(summary$species_taxid==taxid_species[i]))
  }
  row<-summary[ind,]
  
  ind<-which(row$assembly_level=="Complete Genome")
  if(length(ind)>0){
    row<-row[ind,]
  }

  strain_name<-c()
  for(i in 1:nrow(row)){
    strain_name[i]<-paste(row$organism_name[i],gsub("="," ",row$infraspecific_name[i]))
    strain_name[i]<-gsub(" ","_",strain_name[i])
    strain_name[i]<-gsub("\\.","",strain_name[i])
  }
  unique_strain_name<-unique(strain_name)
    
  for(i in 1:length(unique_strain_name)){
    ind<-which(unique_strain_name[i]==strain_name)
    if(length(ind)>1){
      ind<-ind[1]
    }
    ftp_path<-row$ftp_path[ind]
    asm<-row$asm_name[ind]
    download(ftp_path,tmp_genome_path,asm)
  }
}
