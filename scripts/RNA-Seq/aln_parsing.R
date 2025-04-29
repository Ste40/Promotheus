
library(Rsubread)
library(stringr)
library(ProGenome)


root_path="/mnt/e/Stefano40/Setup_env"

setwd(root_path)

Coli_dir <- list.dirs(recursive = FALSE)
Coli_dir <- Coli_dir[ grepl("Escherichia coli", Coli_dir) ]
Coli_dir<- gsub("\\./", "", Coli_dir)
Coli_dir<- paste0(root_path,"/", Coli_dir)

for(i in 1:length(Coli_dir)){
  
  setwd(root_path)
  Coli_subdir <- list.dirs(path=Coli_dir, recursive = TRUE)[-1]
  Coli_subdir<- gsub("\\./", "", Coli_subdir) 
  
  for(j in 1:length(Coli_subdir)){
    
  setwd(Coli_subdir[j])
  
  
  if(identical(list.files(path=root_path, pattern="my_index"), character(0))){
    
    setwd(root_path)
    
    genome <- list.files(path=root_path, pattern= ".fna.gz$", full.names=TRUE, recursive=FALSE)
    buildindex(basename="my_index",reference= genome)
    
  }
    
  
  if(identical(list.files(path=root_path, pattern="*gff.bed$"), character(0))){
    
    setwd(root_path)
    
    GFFbed<- list.files(path=root_path, pattern="*.gff.gz$")
    
    string<- paste('gunzip -k', GFFbed)
    bash<- paste('bash -c', shQuote(string)) 
    system(bash)
    
    GFFbed<- list.files(pattern="*.gff$")
    
    string<- paste('gff2bed <',GFFbed,'> gff.bed')
    bash<- paste('bash -c', shQuote(string))
    system(bash)
    
    file.copy("gff.bed", Coli_subdir[j]) 
    
  } else {
    
    file.copy(paste0(root_path,"/gff.bed"), Coli_subdir[j]) 
  }
  
  
  # aligning reads files fastq to fna genome file
  
  if(identical(list.files(path=Coli_subdir[j], pattern="*fastq.gz.subread.BAM$"), character(0))){
    
    setwd(Coli_subdir[j])
    
    readsfiles= list.files(path=Coli_subdir[j], pattern="fastq.gz$", full.names = FALSE, recursive=FALSE)
    
    Acclist= unique(gsub("\\_.*", "", readsfiles))
    
    if(length(readsfiles)==length(Acclist)*2){
      
      for(k in seq(from=1, to=length(readsfiles), by=2)){
        align(index= file.path(path=root_path, "my_index"), readfile1= readsfiles[k], readfile2=readsfiles[k+1], nthread=THREADS, type="rna", output_file = paste(readsfiles[k],"subread", output_format="BAM" ,sep="."))
        
      }
      
    } else {
      
      for(k in 1:length(Acclist)){  
        align(index= file.path(path=root_path, "my_index"), readfile1= readsfiles[k], nthread=THREADS, type="rna", output_file = paste(readsfiles[k],"subread", output_format="BAM" ,sep="."))
        
      }
    }
    
  }
    
    
    if(identical(list.files(pattern="*sorted.bam$"), character(0))){
      
      
      bamfiles <- list.files(pattern="*BAM$", full.names=FALSE, recursive=FALSE)[1]
      run_name<-gsub("\\..*","",bamfiles)
      
      
      string_run_filtered<- paste0(run_name,'.filtered.bam')
      string_run_sort<- paste0(run_name,'.sorted.bam')
      
      string<- paste('samtools view -@',THREADS,'-q 20 -F4 -o' ,string_run_filtered, bamfiles) 
      bash<- paste('bash -c', shQuote(string))
      system(bash)
      
      string<- paste('samtools sort -@',THREADS,' -o', string_run_sort, string_run_filtered) 
      bash<- paste('bash -c', shQuote(string))
      system(bash)
      
      
      string<- paste('samtools index',string_run_sort)
      bash<- paste('bash -c', shQuote(string))
      system(bash)
      
      file.remove(string_run_filtered)
      
    
    layout<- c()
    
    test<- list.files(pattern="*.sorted.bam$")[1]
    bed<- list.files(path= Coli_subdir[j], pattern="*.bed$", full.names=FALSE, recursive=FALSE)
    string<- paste('infer_experiment.py -s 5000000 -i',test,"sorted.bam -r gff.bed")
    bash<- paste('bash -c', shQuote(string))
    stats<-system(bash, intern=TRUE)
    
    layout[1]<- gsub(".*\\/", "",Coli_subdir[j])
    layout[2]<- stats[grep("This is", stats)]
    
    
    if(layout[2]=="This is PairEnd Data"){
      fraction_stranded<- as.numeric(sub('.*\\:', '', stats[grep("1++", stats, fixed=TRUE)]))			
      fraction_reversely<- as.numeric(sub('.*\\:', '', stats[grep("1+-", stats, fixed=TRUE)]))	
      
    } else {
      fraction_stranded<- as.numeric(sub('.*\\:', '', stats[grep("++", stats, fixed=TRUE)]))			
      fraction_reversely<- as.numeric(sub('.*\\:', '', stats[grep("+-", stats, fixed=TRUE)]))	
      
    }
    
    strand<- fraction_stranded/fraction_reversely
    layout[3]<- strand
    
  
  txt <- paste(layout, collapse = ";")
  txt_path="/mnt/e/Stefano40/Setup_env/layouts.txt"
  write(txt, file = txt_path, append = TRUE)
  
  }
  
  #bam file sorting, filtering and indexing with samtools
  
  if(length(list.files(pattern="*sorted.bam$")) != 
     length(list.files(pattern="*subread.BAM$", full.names=FALSE, recursive=FALSE))){
    
    setwd(Coli_subdir[j])
    bamfiles <- list.files(pattern="*subread.BAM$", full.names=FALSE, recursive=FALSE)[-1]
    
    for (k in 1:length(bamfiles)){
      
      run_name<-gsub("\\..*","",bamfiles[k])
      string_run_filtered<- paste0(run_name,'.filtered.bam')
      string_run_sort<- paste0(run_name,'.sorted.bam')
      
      string<- paste('samtools view -@',THREADS,' -q 20 -F 4 -o' ,string_run_filtered, bamfiles[k]) 
      bash<- paste('bash -c', shQuote(string))
      system(bash)
      
      string<- paste('samtools sort -@',THREADS,' -o', string_run_sort, string_run_filtered) 
      bash<- paste('bash -c', shQuote(string))
      system(bash)
      
      
      string<- paste('samtools index',string_run_sort)
      bash<- paste('bash -c', shQuote(string))
      system(bash)
      
      file.remove(string_run_filtered)
      
    }
    
    if(identical(list.files(path= root_path, pattern="*ptt.ptt$"), character(0))){
      
      setwd(root_path)
      ft <- read.gff(list.files(path= root_path, pattern="_feature_table.txt.gz$", full.names=TRUE, recursive=FALSE))
      
      ptt <- ExtractPtt(ft)
      rnt <- ExtractRnt(ft)
      
      ptt$Length<- round(ptt$Length)
      
      write.ptt(ptt, 'ptt.ptt')
      write.rnt(rnt, 'rnt.rnt')
      
    }
    
  }
  
  }
  
}

