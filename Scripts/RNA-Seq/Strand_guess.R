# This function is necessary to guess reads orientation for TSS files parsing
# Normally this is not necessary, the orientation has been already guessed in file parsing function
# This will be executed in case NDE table has been already generated but something failed in TSS file parsing (eg. cluster execution)

strand_guess <- function(m.o.subpath, GSEdirectory){
  
  layouts<- c()
  
  setwd(GSEdirectory)
  
  test<- list.files(pattern="*.sorted.bam$")[1]
  bed<- list.files(path= GSEdirectory, pattern="*.bed$", full.names=FALSE, recursive=FALSE)
  string<- paste('infer_experiment.py -s 5000000 -i',test,"sorted.bam -r gff.bed")
  bash<- paste('bash -c', shQuote(string))
  stats<-system(bash, intern=TRUE)

  layouts[[1]]<- stats[grep("This is", stats)]
  
  
  if(layouts[[1]]=="This is PairEnd Data"){
    fraction_stranded<- as.numeric(sub('.*\\:', '', stats[grep("1++", stats, fixed=TRUE)]))			
    fraction_reversely<- as.numeric(sub('.*\\:', '', stats[grep("1+-", stats, fixed=TRUE)]))	
    
  } else {
    fraction_stranded<- as.numeric(sub('.*\\:', '', stats[grep("++", stats, fixed=TRUE)]))			
    fraction_reversely<- as.numeric(sub('.*\\:', '', stats[grep("+-", stats, fixed=TRUE)]))	
    
  }
  
  strand<- fraction_stranded/fraction_reversely
  names(layouts[[1]])<- GSE
  layouts[[1]][2]<- strand
  
  strandness<- layouts
  return (strandness)
  rm(layouts)

}