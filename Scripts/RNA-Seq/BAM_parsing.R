# Reject ambiguous layouts before indexing or launching an aligner.
validated_fastq_layout <- function(directory, accessions) {
  if (!length(accessions) || anyNA(accessions) || anyDuplicated(accessions) ||
      any(!grepl("^[SED]RR[0-9]+$", accessions))) stop("Invalid or duplicate run accessions")
  expected <- character()
  layouts <- vapply(accessions, function(code) {
    single <- paste0(code, ".fastq.gz")
    paired <- paste0(code, c("_1.fastq.gz", "_2.fastq.gz"))
    has_single <- file.exists(file.path(directory, single))
    has_pairs <- file.exists(file.path(directory, paired))
    if (has_single && !any(has_pairs)) files <- single
    else if (!has_single && all(has_pairs)) files <- paired
    else stop("Missing, incomplete or ambiguous FASTQ layout for ", code)
    if (any(file.info(file.path(directory, files))$size == 0)) stop("Empty FASTQ for ", code)
    expected <<- c(expected, files)
    length(files)
  }, integer(1))
  if (length(unique(layouts)) != 1L) stop("Mixed single-end and paired-end runs within a GSE are not supported")
  if (!setequal(list.files(directory, pattern = "\\.fastq\\.gz$"), expected))
    stop("Unexpected FASTQ files not associated with requested runs")
  if (layouts[[1]] == 2L) "paired" else "single"
}

run_samtools <- function(arguments, runner = system2) {
  status <- tryCatch(
    suppressWarnings(runner("samtools", args = arguments)),
    error = function(e) 127L
  )
  if (!identical(as.integer(status), 0L)) {
    stop("samtools failed with exit code ", status, ": ", paste(arguments, collapse=" "))
  }
  invisible(TRUE)
}

infer_experiment_stats <- function(bam_file, runner = system2) {
  output <- runner(
    "infer_experiment.py",
    args = c("-s", "5000000", "-i", shQuote(bam_file), "-r", shQuote("gff.bed")),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("infer_experiment.py failed for ", bam_file, " with exit code ", status)
  }
  output
}

BAM_parsing <- function(GSE, m.o.subpath, GSEdirectory, Acclist, threads){

  read_layout <- validated_fastq_layout(GSEdirectory, Acclist)
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

    if(read_layout == "paired"){

    for(i in seq(from=1, to=length(readsfiles), by=2)){
      align(index= file.path(path=m.o.subpath, "my_index"), readfile1= readsfiles[i], readfile2=readsfiles[i+1], nthread=threads, type="rna", output_file = paste(readsfiles[i],"subread", output_format="BAM" ,sep="."))

    }

  } else {

    for(i in seq_along(Acclist)){
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
    run_samtools(c("view", "-@", as.character(threads), "-q", "20", "-o",
                   shQuote(string_run_filter), shQuote(bamfiles)))

    string_run_sort<- paste0(run_name,'.sorted.bam')
    run_samtools(c("sort", "-@", as.character(threads), "-o",
                   shQuote(string_run_sort), shQuote(string_run_filter)))

    run_samtools(c("index", "-@", as.character(threads), shQuote(string_run_sort)))
    unlink(string_run_filter)

  }


    test<- list.files(pattern="*.sorted.bam$")[1]
    stats <- infer_experiment_stats(test)

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

    for (i in seq_along(bamfiles)){

      run_name[i]<-gsub("\\..*","",bamfiles[i])

      string_run_sort<- paste0(run_name[i],'.filtered.bam')
      run_samtools(c("view", "-@", as.character(threads), "-q", "20", "-o",
                     shQuote(string_run_sort), shQuote(bamfiles[i])))

      string_sort<- paste0(run_name[i],'.sorted.bam')
      run_samtools(c("sort", "-@", as.character(threads), "-o",
                     shQuote(string_sort), shQuote(string_run_sort)))

      run_samtools(c("index", "-@", as.character(threads), shQuote(string_sort)))
      unlink(string_run_sort)

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
}
