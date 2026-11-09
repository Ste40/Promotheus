args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else '.')
source(file.path(root,'Scripts/RNA-Seq/MetaData_retrieval.R'))
source(file.path(root,'Scripts/RNA-Seq/NDE_detection.R'))
source(file.path(root,'Scripts/RNA-Seq/BAM_parsing.R'))
fails <- function(expr, pattern) {
  message <- tryCatch({force(expr); ''}, error=function(e) conditionMessage(e))
  stopifnot(grepl(pattern,message,fixed=TRUE))
}
tab <- data.frame(GSM=c('GSM2','GSM1'),SRR=c('SRR2','SRR1'),Group=c('B','A'),GroupNumber=c(2,1))
ordered <- canonical_run_conditions(tab)
scope_metadata <- new.env(parent=globalenv())
sys.source(file.path(root,'Scripts/RNA-Seq/MetaData_retrieval.R'),scope_metadata)
scope_metadata$source <- function(...) invisible(NULL)
scope_metadata$str_split <- strsplit
scope_metadata$Groups_retrieval <- function(gse) tab
retrieved <- scope_metadata$MetaData_retrieval(data.frame(GSE='123',taxon='Toy'),1,data.frame())
stopifnot(retrieved$Accession_codes=='SRR1,SRR2',retrieved$Experimental_condition=='A,B',retrieved$Groups=='1,2')
stopifnot(identical(as.character(ordered$Group),c('A','B')))
dup <- tab;dup$GSM[2] <- dup$GSM[1]
fails(canonical_run_conditions(dup),'one run per GSM')
metadata <- match_count_metadata(c('SRR1.fastq.gz counts.txt','SRR2_1.fastq.gz counts.txt'),
  c('SRR2','SRR1'),c('B','A'),c(2,1))
stopifnot(identical(metadata$conditions,c('A','B')),identical(metadata$groups,c(1,2)))
fails(match_count_metadata('SRR3.fastq.gz counts.txt',c('SRR1'), 'A',1),'do not match')
fails(match_count_metadata(c('SRR1_1.fastq.gz counts.txt','SRR1_2.fastq.gz counts.txt'),c('SRR1','SRR2'),c('A','B'),c(1,2)),'do not match')
run <- function() {
  temp <- tempfile();dir.create(temp);old <- getwd()
  on.exit({setwd(old);unlink(temp,recursive=TRUE)})
  for (file in c('SRR1_1.fastq.gz','SRR1_2.fastq.gz','SRR2.fastq.gz')) writeLines('toy',file.path(temp,file))
  fails(validated_fastq_layout(temp,c('SRR1','SRR2')),'Mixed single-end')
  # The real entrypoint must stop before buildindex/align can run.
  fails(BAM_parsing('GSE1',temp,temp,c('SRR1','SRR2'),1),'Mixed single-end')
  unlink(file.path(temp,'SRR2.fastq.gz'))
  for (file in c('SRR2_1.fastq.gz','SRR2_2.fastq.gz')) writeLines('toy',file.path(temp,file))
  stopifnot(validated_fastq_layout(temp,c('SRR2','SRR1'))=='paired')
  unlink(file.path(temp,'SRR2_2.fastq.gz'))
  fails(validated_fastq_layout(temp,c('SRR1','SRR2')),'incomplete')
  unlink(list.files(temp,full.names=TRUE))
  for (i in 1:2) {
    writeLines('toy',file.path(temp,paste0('SRR',i,'.fastq.gz')))
    file.create(file.path(temp,paste0('SRR',i,'.sorted.bam')))
    write.table(data.frame(GeneID='g',Chr='c',Start=1,End=2,Strand='+',Length=2,Counts=i*11),
      file.path(temp,paste0('SRR',i,'.fastq.gz counts.txt')),sep='\t',quote=FALSE,row.names=FALSE)
  }
  stopifnot(validated_fastq_layout(temp,c('SRR1','SRR2'))=='single')
  # Execute the real count-table loading and metadata matching, substituting only DESeq2.
  scope <- new.env(parent=globalenv())
  sys.source(file.path(root,'Scripts/RNA-Seq/NDE_detection.R'),scope)
  scope$str_split <- strsplit
  scope$DESeqDataSetFromMatrix <- function(countData,colData,design) {
    stopifnot(identical(colnames(countData),c('A','B')),identical(as.numeric(countData[1,]),c(11,22)))
    stop('CAPTURE_OK')
  }
  fails(scope$NDE_detection('GSE1',temp,list(c('This is SingleEnd Data',1)),1,temp,
    data.frame(Experimental_condition='B,A',Groups='2,1'),'toy',file.path(temp,'out'),
    'none',FALSE,.05,c('SRR2','SRR1')),'CAPTURE_OK')
}
run()
cat('PASS: RNA run identity, reordered conditions and layout rejection\n')
