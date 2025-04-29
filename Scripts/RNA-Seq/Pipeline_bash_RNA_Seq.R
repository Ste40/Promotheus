

            #######################################
            # Welcome to Ste40cm RNA-Seq Pipeline # 
            #######################################


INPUT_XLSX <- list.files(
  path = "/input",
  pattern = "RNA_Seq",
  full.names = TRUE
)

setwd("/scripts/RNA-Seq") #Pipeline dir

source("/scripts/RNA-Seq/Installation.R")
source("/scripts/RNA-Seq/Data_retrieval.R")
source("/scripts/RNA-Seq/MetaData_retrieval.R")
source("/scripts/RNA-Seq/BAM_parsing.R")
source("/scripts/RNA-Seq/Operons_detection.R")
source("/scripts/RNA-Seq/NDE_detection.R")
source("/scripts/RNA-Seq/merge_operons_table.R")

# Creating output dir

  if (!file.exists("/output/RNA_Seq")) {
    dir.create(file.path("/output/RNA_Seq"), recursive = TRUE)
  }

# Check for package installation

message("Checking if the required packages are installed ")

if(Packages_installation(packages)==FALSE) {
  stop("Oooops! Something went wrong during package installation")
} else {message ("Installation successful!")}

# Loading packages

message("Loading the required libraries..")
suppressWarnings(suppressPackageStartupMessages({
library(Rsubread)
library(GEOquery)
library(biomartr)
library(dplyr)
library(DESeq2)
library(reshape2)
library(stringr)
library(ProGenome)
library(EnvStats)
library(openxlsx)
library(curl)
library(data.table)
library(rtracklayer)
library(ggplot2)
library(Biostrings)
}))

message("Done!")
#######################################
#######################################

microorganism_table= read.xlsx(INPUT_XLSX)

xlsx <- data.frame(m.o.name= character(),
                   GSE= character(), 
                   Experimental_condition= character(), 
                   Groups= character(),
                   Accession_codes= character())


if(!file.exists(file.path("/output/RNA_Seq/Experimental_Groups.xlsx"))){
  message("I'm retrieving metadata, please wait")
  for(g in 1:nrow(microorganism_table)){
    xlsx_app<- suppressWarnings(MetaData_retrieval(microorganism_table, g, xlsx))
    xlsx<- rbind(xlsx, xlsx_app)
  }

  write.xlsx(xlsx, file.path("/output/RNA_Seq/Experimental_Groups.xlsx"), overwrite=TRUE)
  if(nrow(xlsx) != length(str_split(microorganism_table$GSE, ";")[])) {
    message("Oooops! Something went wrong during metadata retrieval, please check 'Experimental_Groups.xlsx' file") 
  } else {
    message ("Success: metadata retrieved!")
  }
} else {
    message ("Success: Experimental Groups already defined!")
}


cat(
'Before proceeding, check if the groups and replices have been assigned correctly within table "Experimental_Groups.xlsx". 
If not, please edit the file and overwrite it. When you are done, press enter:')
invisible(scan("stdin", character(), nlines = 1, quiet = TRUE))


#######################################
#######################################


Conditions<- read.xlsx(list.files(path="/output/RNA_Seq", pattern="*Groups.xlsx", full.names = TRUE)[1])
GSEs <- Conditions[,2] 
gse_count <- aggregate(GSE ~ m.o.name, data = Conditions, FUN = function(x) length(unique(x)))
original_threshold <- operon_threshold
for (r in 1:length(GSEs)){

  temp_threshold <- original_threshold
  GSE<- GSEs[r] 
  Name_ID<- Conditions$m.o.name[r]
  Assembly_v <- microorganism_table$reference_genome[r]
  m.o.correct= str_replace_all(Name_ID, "[[:punct:]]", " ")
  GSEconditions<- Conditions[r,c("Experimental_condition", "Groups")]
  m.o.subpath= file.path("/output/RNA_Seq/",m.o.correct)
  Acclist<- str_split(Conditions$Accession_codes[r], ",")[[1]]
  GSEdirectory<- file.path(m.o.subpath,paste(GSE, m.o.correct))
  threads<- as.numeric(THREADS)
  memory<- as.numeric(MEMORY)
  od_full_tabs= file.path("/output/RNA_Seq/results")
    
  if (!file.exists(file.path(od_full_tabs, paste(GSE,m.o.correct,"full_tab.txt")))){ 
    
    suppressWarnings(dir.create(m.o.subpath))
    
    message(paste("I'm downloading the required files for", Name_ID, GSE))
    retr = Data_retrieval(GSE, Acclist, GSEdirectory, m.o.subpath, args, threads, GSEconditions, Assembly_v)
    if (retr==TRUE) {message ("Success: files downloaded!")}
    if (retr == "skip"){
      next
    }
    
    message(paste("I'm aligning the reads for", GSE, "to reference genome, please wait"))
    strandness<- BAM_parsing(GSEs, r, m.o.subpath, GSEdirectory, Acclist, threads)
    if(strandness[[1]][1] != "This is SingleEnd Data" & strandness[[1]][1] != "This is PairEnd Data") {
      stop("Oooops! Something went wrong during alignment or .bam parsing")
    } else {message ("Success: .BAM file analyzed!")}

    if(OPERONS_DETECTION==TRUE){
      message("I'm detecting operons, please wait")
      if(Operons_detection(m.o.subpath, threads, memory, Acclist, GSEdirectory, GSEconditions, temp_threshold, gse_count) != TRUE) {
        warning(paste("Oooops! Something went wrong during Operon detection for", basename(GSEdirectory), 
        "skipping operons detection for this GSE"))
      } else {message ("Success: Operons detected!")}
    } else {
      message("OPERONS_DETECTION is set to FALSE, operons detection will be skipped")
    }

    message("I'm detecting Not-Differentially Expressed (NDE) genes, please wait")
    if(NDE_detection(GSEdirectory, strandness, threads, m.o.subpath, GSEconditions, m.o.correct, od_full_tabs, CORRECTION_DE, OPERONS_DETECTION, PVALUE) != TRUE) {
      stop("Oooops! Something went wrong during NDE identification")
    } else {message ("Success: NDE genes detected!")}

    if(OPERONS_DETECTION==TRUE){
    tryCatch({
      message("I'm updating the GSE table based on operons information, please wait")
      merge_operons(m.o.subpath, od_full_tabs, GSEdirectory, GSE, m.o.correct)
      message("Success: NDE genes detected!")
    }, error = function(e) {
      stop("Oooops! Something went wrong during operons merge: ", e$message)
    })
    }

    message(paste(GSE, "analisys completed!"))
  
    #unlink(GSEdirectory, force=TRUE, recursive = TRUE)
    
    if(r == length(GSEs)){
      message("Analisys completed! All the GSEs have been analyzed")
    }
  } else {message(paste(GSE,"has been already been analyzed, I proceed with the next one.."))}
}

if (all(sapply(GSEs, function(GSE) file.exists(file.path(od_full_tabs, paste(GSE, m.o.correct, "full_tab.txt")))))) {
  message(paste("All the GSEs have been analyzed!"))
}



