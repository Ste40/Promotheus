
###### Script for fastq download #######

## Nel nostro caso, la tabella si chiama 'Experimental_Groups.xlsx', 
##che si trova nel percorso '/mnt/e/Stefano40/SQ_26_10/Pipeline_RNA-Seq/Experimental_Groups.xlsx' !!!!DA NON CANCELLARE O MODIFICARE !!!!

## la tabella va splittata per ogni PC che useremo, quindi se ne usiamo 2, 
##bisogna creare due tabelle nuove splittando l'originale.

## Di seguito un esempio degli argomenti input dello script 'Rscript download.R <arg1> <arg2> <arg3>' 

#args = c("/mnt/e/Stefano40/Setup_env",     ####Percorso alla cartella root dove salvare i file 
#       "/mnt/e/Stefano40/Setup_env/Experimental_Groups.xlsx",       ####Percorso alla tabella, ricordare l'estensione .xlsx
#    "4")                                  ####numero di threads per il download, controllate il numero di cores per CPU, poi fate ncores*2


list.of.packages = c("openxlsx", "stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(openxlsx)
library(stringr)

args = commandArgs(trailingOnly=TRUE) 

Conditions<- read.xlsx(list.files(path=args[1], pattern="*Groups.xlsx", full.names = TRUE))
GSEs <- Conditions[,2] 

for (r in 1:length(GSEs)){
  
  GSE<- GSEs[r] 
  Name_ID<- Conditions$m.o.name[r]
  m.o.correct= str_replace_all(Name_ID, "[[:punct:]]", " ")
  GSEconditions<- Conditions[r,c("Experimental_condition", "Groups")]
  m.o.subpath= file.path(args[1],m.o.correct)
  threads<- as.numeric(args[3])
  
  Acclist<- str_split(Conditions$Accession_codes[r], ",")[[1]]
  GSEdirectory<- file.path(m.o.subpath,paste(GSE, m.o.correct))
    
  suppressWarnings(dir.create(m.o.subpath))
    
  message(paste("I'm downloading the required files for", Name_ID, GSE))
  
  token<- c()
  
  #if(length(str_split(GSEconditions$Groups,",")[[1]])==length(Acclist)){
    
    suppressWarnings(dir.create(GSEdirectory))  
    setwd(GSEdirectory)
    
    if(identical(list.files(path = GSEdirectory, pattern="*.fastq.gz$"), character(0))){
      
      message("I'm downloading reads Fastq files, please wait. This operation may require time depending on internet connection and NCBI server status..")
      
      for (i in 1:length(Acclist)){
        
        message(paste("I'm downloading Fastq file/s for", Acclist[i]))
        #cmd = paste('fastq-dump --gzip --split-files --readids --skip-technical --dumpbase --clip -O',shQuote(GSEdirectory), Acclist[i]) @not parallelized
        cmd = paste('parallel-fastq-dump -s', Acclist[i] ,'-t', threads, '--gzip --skip-technical --split-files --readids --dumpbase --clip --tmpdir',shQuote(GSEdirectory)) #parallelized
        system(cmd)
        
      }
    } 
    
    #if(length(list.files(path = GSEdirectory, pattern="*.fastq.gz$")) != length(Acclist) | 
     #  length(list.files(path = GSEdirectory, pattern="*.fastq.gz$")) != length(Acclist)*2){
      
     # keep<- list.files(path = GSEdirectory, pattern="*_1.fastq.gz|*_2.fastq.gz", full.names = TRUE)
      #junk<- list.files(GSEdirectory, "*.fastq.gz", full.names = TRUE)
      
      #'%ni%' <- Negate('%in%')
      #junk<- junk[junk %ni% keep]
      #invisible(file.remove(junk))
      
    #}
    
    #if(length(list.files(path = GSEdirectory, pattern="*.fastq.gz$")) == length(Acclist) | 
     #  length(list.files(path = GSEdirectory, pattern="*.fastq.gz$")) == length(Acclist)*2) {
      
      token[1]<- TRUE
      
   # } else {
      
    #  token[1]<- FALSE
      
    #}
    
    setwd(m.o.subpath)
  }# else {token[1]<- FALSE}
  
    if (token==TRUE) {message ("Success: files downloaded!")}
    if (token == FALSE){
      message(paste("Something went wrong during fastq download for", GSE, "I proceed with the next one"))
      next

  }
}