


BAM_parsing = function(GSEs, r, m.o.subpath, GSEdirectory, Acclist, threads){
  
  setwd(m.o.subpath)
  layouts<- list()
  
# First task: alignment
# building index for genome
  
  if(identical(list.files(path=m.o.subpath, pattern="my_index"), character(0))){

    genome <- list.files(path=m.o.subpath, pattern= ".fna.gz$", full.names=TRUE, recursive=FALSE)
    buildindex(basename="my_index",reference= genome)
    
  }
   
# coverting GFF to bed (input for infer_experiment.py)

  if(identical(list.files(pattern="*gff.bed$"), character(0))){
    
    GFFbed<- list.files(pattern="*.gff.gz$")

    string<- paste('gunzip -k', GFFbed)
    bash<- paste('bash -c', shQuote(string)) 
    system(bash)
    
    GFFbed<- list.files(pattern="*.gff$")

    string<- paste('gff2bed <',GFFbed,'> gff.bed')
    bash<- paste('bash -c', shQuote(string))
    system(bash)
    
    file.copy("gff.bed", GSEdirectory) 
    
  } else {
    
    file.copy("gff.bed", GSEdirectory) 
  }
  
  setwd(GSEdirectory)
  
# aligning reads files fastq to fna genome file

  if(identical(list.files(pattern="*fastq.gz.subread.BAM$"), character(0))){
    readsfiles <- list.files(pattern="*.fastq.gz$", full.names=FALSE, recursive=FALSE)
  
    if(length(readsfiles)==length(Acclist)*2){
    
    for(i in seq(from=1, to=length(readsfiles), by=2)){
      align(index= file.path(path=m.o.subpath, "my_index"), readfile1= readsfiles[i], readfile2=readsfiles[i+1], nthread=threads, type="rna", output_file = paste(readsfiles[i],"subread", output_format="BAM" ,sep="."))
      
    }
      
  } else {
    
    for(i in 1:length(Acclist)){  
      align(index= file.path(path=m.o.subpath, "my_index"), readfile1= readsfiles[i], nthread=threads, type="rna", output_file = paste(readsfiles[i],"subread", output_format="BAM" ,sep="."))
      
      }
    }
  }
  
########################################
########################################
  
  
#checking if the library is fr or rf with infer_experiment.py and samtools

  if(identical(list.files(pattern="*sorted.bam$"), character(0))){
    
    bamfiles <- list.files(pattern="*BAM$", full.names=FALSE, recursive=FALSE)[1]
    run_name<-gsub("\\..*","",bamfiles)
    
    string_run_filter<- paste0(run_name,'.filtered.bam')
    string<- paste('samtools view -q 20 -o',string_run_filter, bamfiles)
    bash<- paste('bash -c', shQuote(string))
    system(bash)

    string_run_sort<- paste0(run_name,'.sorted.bam')
    string<- paste('samtools sort -o',string_run_sort, string_run_filter)
    bash<- paste('bash -c', shQuote(string))
    system(bash)

    string<- paste('samtools index',string_run_sort)
    bash<- paste('bash -c', shQuote(string))
    system(bash)

  }
  
  
    test<- list.files(pattern="*.sorted.bam$")[1]
    bed<- list.files(path= m.o.subpath, pattern="*.bed$", full.names=FALSE, recursive=FALSE)
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

  

#bam file sorting, filtering and indexing with samtools

  if(length(list.files(path= GSEdirectory, pattern="*sorted.bam$")) != 
            length(list.files(pattern="*subread.BAM$", full.names=FALSE, recursive=FALSE))){
  
    run_name<-c()
    bamfiles <- list.files(pattern="*subread.BAM$", full.names=FALSE, recursive=FALSE)[-1]
  
    for (i in 1:length(bamfiles)){
    
      run_name[i]<-gsub("\\..*","",bamfiles[i])
    
      string_run_sort<- paste0(run_name[i],'.filtered.bam')
      string<- paste('samtools view -q 20 -o',string_run_sort, bamfiles[i])
      bash<- paste('bash -c', shQuote(string))
      system(bash)
    
      string_sort<- paste0(run_name[i],'.sorted.bam')
      string<- paste('samtools sort -o',string_sort, string_run_sort)
      bash<- paste('bash -c', shQuote(string))
      system(bash)
    
      string<- paste('samtools index',string_sort)
      bash<- paste('bash -c', shQuote(string))
      system(bash)
    
    }
  }
  
#creating ptt and rnt (rockhopper's input files) from feature table
  
  if(identical(list.files(path= m.o.subpath, pattern="*ptt.ptt$"), character(0))){
    
    setwd(m.o.subpath)
    ft <- read.gff(list.files(path= m.o.subpath, pattern="_feature_table.txt.gz$", full.names=TRUE, recursive=FALSE))
    
    ptt <- ExtractPtt(ft)
    rnt <- ExtractRnt(ft)
    
    ptt$Length<- round(ptt$Length)
    
    write.ptt(ptt, 'ptt.ptt')
    write.rnt(rnt, 'rnt.rnt')
    
  }
  
  
  strandness<- layouts
  return (strandness)
  rm(layouts)
}
