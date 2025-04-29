
#This function downloads the necessary files for data processing 

Annotation_retrieval <- function(ftp_url){

  ftp_path<- ftp_url
  ftp_string<- gsub(".*/", "", ftp_path)
  genome_ftp<- paste0(ftp_path, "/", ftp_string, "_genomic.fna.gz")
  GFF_ftp<- paste0(ftp_path, "/", ftp_string, "_genomic.gff.gz")
  ft_ftp<- paste0(ftp_path, "/", ftp_string, "_feature_table.txt.gz")
  
  genome<- curl_download(genome_ftp, destfile=paste0(m.o.subpath,"/", ftp_string, "_genomic.fna.gz"), mode="wb")
  GFF<- curl_download(GFF_ftp, destfile=paste0(m.o.subpath,"/", ftp_string, "_genomic.gff.gz"), mode="wb")
  ft<- curl_download(ft_ftp, destfile=paste0(m.o.subpath,"/", ftp_string, "_feature_table.txt.gz"), mode="wb")

} 



Data_retrieval <- function(GSE, Acclist, GSEdirectory, m.o.subpath, args, threads, GSEconditions, Assembly_v){
  
  if(length(str_split(GSEconditions$Groups,",")[[1]])==length(Acclist)){

    token<- c()
    suppressWarnings(dir.create(GSEdirectory))  
    setwd(GSEdirectory)
    
    fq_files_check <- sapply(Acclist, function(code) {
    any(file.exists(list.files(path = GSEdirectory, pattern = paste0(code, "_[12]?\\.fastq\\.gz$"))))
    })

    if(!any(fq_files_check)){
      
      message("I'm downloading reads Fastq files, please wait. This operation may require time depending on internet connection and NCBI server status..")
      
      for (i in 1:length(Acclist)){
          
        if(all(!file.exists(list.files(path = GSEdirectory, pattern = paste0(Acclist[i], "_[12]?\\.fastq\\.gz$|", Acclist[i], "\\.fastq\\.gz$"))))){

          message(paste("I'm downloading Fastq file/s for", Acclist[i]))
          #cmd = paste('fastq-dump --gzip --split-files --readids --skip-technical --dumpbase --clip -O',shQuote(GSEdirectory), Acclist[i]) @not parallelized
          cmd = paste('fasterq-dump -m 1G -p --split-3 --skip-technical -x', Acclist[i] ,'-e', threads, '-t',shQuote(GSEdirectory)) #parallelized
          system(cmd)

          if (length(list.files(path = GSEdirectory, pattern=paste0(Acclist[i],".*\\.fastq$")))>1){
            message(paste(Acclist[i],"contains paired-end reads!"))

            keep<- list.files(path = GSEdirectory, pattern=paste0(Acclist[i],"*_1.fastq|*_2.fastq|*.fastq.gz$"), full.names = TRUE)
            junk<- list.files(GSEdirectory, full.names = TRUE)
        
            '%ni%' <- Negate('%in%')
            junk<- junk[junk %ni% keep]
            invisible(file.remove(junk))
    
            message(paste(Acclist[i], "files donwloaded! I'm compressing them to save space.."))
            to_compress <- list.files(path = GSEdirectory, pattern=paste0(Acclist[i],".*\\.fastq$"))
            for (k in 1:length(to_compress)){
              cmd <- paste('gzip', to_compress[k])
              system(cmd) 
            }
          } else if(length(list.files(path = GSEdirectory, pattern=paste0(Acclist[i],".*\\.fastq$")))==1){
            message(paste(Acclist[i],"contains single-end reads!"))

            keep<- list.files(path = GSEdirectory, pattern=paste0(Acclist[i],"*.fastq|*.fastq.gz$"), full.names = TRUE)
            junk<- list.files(GSEdirectory, full.names = TRUE)
        
            '%ni%' <- Negate('%in%')
            junk<- junk[junk %ni% keep]
            invisible(file.remove(junk))
    
            message(paste(Acclist[i], "files donwloaded! I'm compressing them to save space.."))
            to_compress <- list.files(path = GSEdirectory, pattern=paste0(Acclist[i],".*\\.fastq$"))
            for (k in 1:length(to_compress)){
              cmd <- paste('gzip', to_compress[k])
              system(cmd) 
            }
          } else {
            warning(paste("Download of", Acclist[i], "failed"))
          }
        } else {
           message(paste(Acclist[i], "already downloaded!"))
        }
      }
    } 
   
    # if(length(list.files(path = GSEdirectory, pattern="*.fastq.gz$")) != length(Acclist) | 
    #     length(list.files(path = GSEdirectory, pattern="*.fastq.gz$")) != length(Acclist)*2){
        
    #       keep<- list.files(path = GSEdirectory, pattern="*_1.fastq.gz|*_2.fastq.gz", full.names = TRUE)
    #       junk<- list.files(GSEdirectory, full.names = TRUE)
        
    #       '%ni%' <- Negate('%in%')
    #       junk<- junk[junk %ni% keep]
    #       invisible(file.remove(junk))
  
    #     }

    if(length(list.files(path = GSEdirectory, pattern="*.fastq.gz$")) == length(Acclist) | 
           length(list.files(path = GSEdirectory, pattern="*.fastq.gz$")) == length(Acclist)*2) {
        
          token[1]<- TRUE
        
      } else {
      
        token[1]<- FALSE
      
      }

    setwd(m.o.subpath)
    
    
  # downloading the genome fna file, GFF and feature table, if don't exist in m.o root directory.
    
    message("I'm searching genome files on databases, please wait.") 
    
    if(identical(list.files(pattern="*.fna.gz$|*.gff.gz$|*_feature_table.txt.gz$"), character(0))){
      if (identical(list.files(path="/output/",pattern="*assembly_summary.txt$"), character(0))){
        
        message("Downloading assembly summary from RefSeq..")
        summary<- curl_download("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"
                                , destfile=file.path("/output/", "assembly_summary.txt"), mode="wb")
      }
        
      summary<- fread(file.path("/output/","assembly_summary.txt"), quote=F)
      ftp_url<- summary[summary$asm_name==Assembly_v]$ftp_path
      
      if(identical(ftp_url, character(0))){
        
        message("I can't find any ftp link for the provided genome assembly accession code, I'll search for organism name../n")
        Nameparse<- strsplit(Name_ID, " ")[[1]]
        Nameparse_1<- paste(Nameparse[1],Nameparse[2])
        Nameparse_2<- tolower(Nameparse[length(Nameparse)])
        
        summary_red<-summary[agrep(Nameparse_1, summary$organism_name, fixed=F, ignore.case=TRUE),]
        summary_red_ref <- summary_red[summary_red$refseq_category=="reference genome",]

        if(nrow(summary_red_ref) != 0){
          if(nrow(summary_red_ref) > 0){
            ftp_url<- summary_red_ref$ftp_path[summary_red_ref$organism_name==Name_ID]
          } else {
          ftp_url<- summary_red_ref$ftp_path
          }
        } else {
          summary_red_ref <- summary_red[summary_red$refseq_category=="representative genome",]
          ftp_url<- summary_red_ref$ftp_path
        }

        if(identical(ftp_url, character(0))){

          summary_red$infraspecific_name <- gsub("strain=","",summary_red$infraspecific_name)
          summary_red_spec<-summary_red[summary_red$infraspecific_name %like% Nameparse_2, ]
          summary_red_spec_genome<- summary_red_spec[summary_red_spec$refseq_category=="reference genome", ]

        if(is.na(summary_red_spec_genome[1,1])){
        
          ftp_url<-summary_red_spec[summary_red_spec$assembly_level=="Complete Genome",][1]$ftp_path
           
        } else {
        
          ftp_url<-summary_red_spec_genome[summary_red_spec_genome$assembly_level=="Complete Genome",][1]$ftp_path
          
        } 
      }
      if(identical(ftp_url, character(0))){
        
        message(paste("An error occurred during annotation retrieval, please check m.o.name for", GSE))
        token[2]=FALSE
        
        }
      }
      
      annotations= Annotation_retrieval(ftp_url)
      
      
      if(identical(list.files(pattern="*.fna.gz$|*.gff.gz$|*_feature_table.txt.gz$"), character(0))){
        
        token[2]=FALSE
        
      }
      
      token[2]<- TRUE
     
    } else {
      
      message(paste("Success: annotation files already locally stored!"))
      genome<- list.files(pattern="*.fna.gz$")
      token[2]<- TRUE
      
    }
    


  if(all(token)){
    
    message=TRUE
    
  } else { 
  
    message(paste("The files are incomplete or damaged", GSE, "will be skipped"))
    message="skip"
    unlink(GSEdirectory, recursive=TRUE)
    
  }
    
  } else {
    
    message(paste("The number of experimental conditions is different from the number of runs, the analysis of the experiment", GSE, "will be skipped"))
    message="skip"
  }
    
  return(message)
  
}

