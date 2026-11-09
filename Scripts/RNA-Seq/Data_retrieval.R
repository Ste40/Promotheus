
# This function downloads the necessary files for data processing.

refseq_file_location <- function(ftp_path, suffix) {
  if (length(ftp_path) != 1 || is.na(ftp_path) || !nzchar(ftp_path) || ftp_path == "na") {
    stop("RefSeq assembly has no usable ftp_path")
  }

  ftp_path <- sub("/+$", "", ftp_path)
  ftp_string <- basename(ftp_path)
  if (!nzchar(ftp_string)) {
    stop("Unable to determine RefSeq assembly name from ftp_path: ", ftp_path)
  }

  list(
    url = paste0(ftp_path, "/", ftp_string, suffix),
    filename = paste0(ftp_string, suffix)
  )
}

Annotation_retrieval <- function(ftp_url, destination_dir) {

  genome_location <- refseq_file_location(ftp_url, "_genomic.fna.gz")
  gff_location <- refseq_file_location(ftp_url, "_genomic.gff.gz")
  feature_location <- refseq_file_location(ftp_url, "_feature_table.txt.gz")

  genome <- curl_download(
    genome_location$url,
    destfile = file.path(destination_dir, genome_location$filename),
    mode = "wb"
  )
  GFF <- curl_download(
    gff_location$url,
    destfile = file.path(destination_dir, gff_location$filename),
    mode = "wb"
  )
  ft <- curl_download(
    feature_location$url,
    destfile = file.path(destination_dir, feature_location$filename),
    mode = "wb"
  )

}

fastq_files_for_accession <- function(accession, directory, compressed = FALSE) {
  extension <- if (compressed) "\\.fastq\\.gz$" else "\\.fastq$"
  pattern <- paste0("^", accession, "(_[123])?", extension)
  list.files(path = directory, pattern = pattern, full.names = TRUE)
}

normalize_split_fastq_files <- function(accession, directory) {
  files <- fastq_files_for_accession(accession, directory, compressed = FALSE)
  pair_paths <- file.path(directory, paste0(accession, c("_1.fastq", "_2.fastq")))
  if (all(file.exists(pair_paths))) {
    auxiliary <- setdiff(files, pair_paths)
    if (length(auxiliary) > 0L) {
      unlink(auxiliary, force = TRUE)
    }
    return(pair_paths)
  }
  files
}

fastq_accession_complete <- function(accession, directory, compressed = TRUE) {
  files <- fastq_files_for_accession(accession, directory, compressed = compressed)
  length(files) %in% c(1L, 2L) && all(file.info(files)$size > 0)
}

compress_fastq_files <- function(files, threads, runner = system2) {
  if (length(files) == 0) {
    return(FALSE)
  }

  compress_one <- function(path) {
    tryCatch(
      as.integer(runner("gzip", args = c("-1", "-f", "--", shQuote(path)))),
      error = function(e) 127L
    )
  }
  workers <- max(1L, min(as.integer(threads), length(files)))
  if (workers > 1L && .Platform$OS.type != "windows") {
    statuses <- unlist(parallel::mclapply(files, compress_one, mc.cores=workers))
  } else {
    statuses <- vapply(files, compress_one, integer(1))
  }
  all(statuses == 0L)
}

download_fastq_with_retry <- function(accession, directory, threads,
                                      max_attempts = 3L, base_delay = 30,
                                      runner = system2, sleep_fn = Sys.sleep) {
  if (max_attempts < 1) {
    stop("max_attempts must be at least 1")
  }

  for (attempt in seq_len(max_attempts)) {
    partial_files <- c(
      fastq_files_for_accession(accession, directory, compressed = FALSE),
      fastq_files_for_accession(accession, directory, compressed = TRUE)
    )
    if (length(partial_files) > 0) {
      unlink(partial_files, force=TRUE)
    }

    message(paste("I'm downloading Fastq file/s for", accession,
                  paste0("(attempt ", attempt, "/", max_attempts, ")")))
    arguments <- c(
      "-m", "1G", "-p", "--split-3", "--skip-technical", "-x", accession,
      "-e", as.character(threads), "-t", shQuote(directory)
    )
    status <- tryCatch(
      suppressWarnings(runner("fasterq-dump", args = arguments)),
      error = function(e) {
        warning(paste("Unable to start fasterq-dump for", accession, ":", e$message))
        127L
      }
    )
    if (identical(as.integer(status), 0L)) {
      normalize_split_fastq_files(accession, directory)
    }
    downloaded_files <- fastq_files_for_accession(accession, directory, compressed = FALSE)
    valid_files <- fastq_accession_complete(accession, directory, compressed = FALSE)
    if (identical(as.integer(status), 0L) && valid_files) {
      return(TRUE)
    }

    if (length(downloaded_files) > 0) {
      unlink(downloaded_files, force=TRUE)
    }
    warning(paste("fasterq-dump failed for", accession,
                  "with exit code", status, "on attempt", attempt))
    if (attempt < max_attempts && base_delay > 0) {
      delay <- base_delay * 2^(attempt - 1)
      message(paste("Retrying", accession, "in", delay, "seconds"))
      sleep_fn(delay)
    }
  }
  FALSE
}



Data_retrieval <- function(GSE, Acclist, GSEdirectory, m.o.subpath, args, threads, GSEconditions, Assembly_v,
  Name_ID) {

  if(length(str_split(GSEconditions$Groups,",")[[1]])==length(Acclist)){

    token<- c()
    suppressWarnings(dir.create(GSEdirectory))
    setwd(GSEdirectory)

    fq_files_check <- vapply(Acclist, function(code) {
      fastq_accession_complete(code, GSEdirectory, compressed=TRUE)
    }, logical(1))

    if(!all(fq_files_check)){

      message("I'm downloading reads Fastq files, please wait. This operation may require time depending on internet connection and NCBI server status..")

      missing_accessions <- Acclist[!fq_files_check]
      requested_workers <- suppressWarnings(as.integer(Sys.getenv("MAX_GSE_PARALLEL", "1")))
      if (is.na(requested_workers) || requested_workers < 1L) {
        requested_workers <- 1L
      }
      download_workers <- min(requested_workers, length(missing_accessions), as.integer(threads))
      threads_per_download <- max(1L, floor(as.integer(threads) / download_workers))

      download_one <- function(accession) {
        download_ok <- download_fastq_with_retry(
          accession, GSEdirectory, threads_per_download,
          max_attempts = 3L, base_delay = 30
        )
        if (!download_ok) {
          warning(paste("Download of", accession, "failed after all retry attempts"))
          return(FALSE)
        }

        to_compress <- fastq_files_for_accession(accession, GSEdirectory, compressed = FALSE)
        read_type <- if (length(to_compress) == 2L) "paired-end" else "single-end"
        message(paste(accession, "contains", read_type, "reads!"))
        message(paste(accession, "files downloaded! I'm compressing them to save space.."))
        if (!compress_fastq_files(to_compress, threads_per_download)) {
          warning(paste("gzip failed for one or more FASTQ files of", accession))
          return(FALSE)
        }
        fastq_accession_complete(accession, GSEdirectory, compressed = TRUE)
      }

      if (download_workers > 1L && .Platform$OS.type != "windows") {
        message("Downloading ", length(missing_accessions), " accessions with ", download_workers,
          " workers (", threads_per_download, " threads each)")
        parallel::mclapply(
          missing_accessions,
          download_one,
          mc.cores = download_workers,
          mc.preschedule = FALSE
        )
      } else {
        lapply(missing_accessions, download_one)
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

    complete_accessions <- vapply(Acclist, function(code) {
      fastq_accession_complete(code, GSEdirectory, compressed=TRUE)
    }, logical(1))

    if(all(complete_accessions)) {

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

      annotations <- Annotation_retrieval(ftp_url, m.o.subpath)


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
