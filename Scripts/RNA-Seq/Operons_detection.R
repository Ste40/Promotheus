# Rockhopper uses commas within a replicate group and spaces between groups.
rockhopper_read_groups <- function(directory, accessions, group_ids) {
  accessions <- as.character(accessions)
  group_ids <- as.character(group_ids)
  if (!length(accessions) || length(accessions) != length(group_ids) ||
      anyNA(accessions) || anyNA(group_ids) || anyDuplicated(accessions) ||
      any(!grepl("^[SED]RR[0-9]+$", accessions)) || any(!nzchar(group_ids)))
    stop("Invalid or unmatched run/group metadata for Rockhopper")
  layouts <- integer(length(accessions))
  expected <- character()
  reads <- vapply(seq_along(accessions), function(i) {
    code <- accessions[[i]]
    single <- paste0(code, ".fastq.gz")
    pair <- paste0(code, c("_1.fastq.gz", "_2.fastq.gz"))
    has_single <- file.exists(file.path(directory, single))
    has_pair <- file.exists(file.path(directory, pair))
    if (has_single && !any(has_pair)) files <- single
    else if (!has_single && all(has_pair)) files <- pair
    else stop("Missing or ambiguous FASTQ layout for ", code)
    if (any(file.info(file.path(directory, files))$size == 0)) stop("Empty FASTQ for ", code)
    layouts[[i]] <<- length(files)
    expected <<- c(expected, files)
    paste(files, collapse = "%")
  }, character(1))
  if (length(unique(layouts)) != 1L) stop("Mixed FASTQ layouts are not supported by this adapter")
  if (!setequal(list.files(directory, pattern = "\\.fastq\\.gz$"), expected))
    stop("Unexpected FASTQ files for Rockhopper")
  vapply(unique(group_ids), function(group) paste(reads[group_ids == group], collapse = ","),
         character(1), USE.NAMES = FALSE)
}

Operons_detection <- function(m.o.subpath, threads, memory, Acclist, GSEdirectory, GSEconditions, temp_threshold,
  gse_count, GSE, m.o.correct) {

  # Estrae il nome del microrganismo dalla directory m.o.subpath
  microorganism <- basename(m.o.subpath)

  # Ottiene il numero di GSE per il microrganismo dal data.frame gse_count
  gse_for_mo <- as.numeric(gse_count$GSE[str_replace_all(gse_count$m.o.name, "[[:punct:]]", " ") == microorganism])
  if(length(gse_for_mo) == 0) {
    warning("Microorganism ", microorganism, " not found in gse_count. Using user-provided operon_threshold.")
  } else if(temp_threshold > gse_for_mo) {
    warning("User-defined operon_threshold (", temp_threshold, ") is greater than available GSE (", gse_for_mo, ") for microorganism ", microorganism, ". Reassigning operon_threshold to ", gse_for_mo, ".")
    temp_threshold <- gse_for_mo
  }

    # Se il file degli operoni non esiste ancora nella directory corrente, procedi con la rilevazione
  if(identical(list.files(path = GSEdirectory, pattern = "*_operons.txt$"), character(0))) {

    Rockhopperpath <- "/scripts/RNA-Seq/Rockhopper.jar"

    if (!is.data.frame(GSEconditions) || nrow(GSEconditions) != 1L ||
        !"Groups" %in% names(GSEconditions)) stop("Missing Rockhopper group assignments")
    group_ids <- strsplit(as.character(GSEconditions$Groups[[1]]), ",", fixed = TRUE)[[1]]
    groupsvector <- rockhopper_read_groups(GSEdirectory, Acclist, group_ids)
    setwd(GSEdirectory)
    if (!length(list.files(m.o.subpath, pattern = "\\.fna$"))) {
      genomes <- list.files(m.o.subpath, pattern = "\\.fna\\.gz$", full.names = TRUE)
      if (!length(genomes)) stop("No reference FASTA for Rockhopper")
      status <- system2("gunzip", c("-k", shQuote(genomes)))
      if (status != 0L) stop("Reference decompression failed for Rockhopper")
    }
    message("Rockhopper replicate groups: ", paste(groupsvector, collapse = " | "))
    status <- system2("java", c(paste0("-Xmx", memory, "m"), "-cp", shQuote(Rockhopperpath),
      "Rockhopper", "-o", shQuote(GSEdirectory), "-g", shQuote(m.o.subpath), shQuote(groupsvector)))
    if (status != 0L) stop("Rockhopper failed with exit code ", status, " for ", GSE)
    if (!length(list.files(GSEdirectory, pattern = "_operons.txt$")))
      stop("Rockhopper completed without an operon table for ", GSE)
  }

  parent_dir <- dirname(GSEdirectory)
  gse_dirs <- list.dirs(path = parent_dir, full.names = TRUE, recursive = FALSE)
  detected_count <- 0
  for(gse in gse_dirs) {
    if(length(list.files(path = gse, pattern = "*_operons.txt$")) > 0) {
      detected_count <- detected_count + 1
    }
  }

mode_filename <- file.path(parent_dir, paste0(microorganism, "_operons_mode.txt"))
if(file.exists(mode_filename)) {
  cat("Operon detection already performed for", temp_threshold,
      "different GSE datasets for this microorganism. Skipping detection.\n")
} else if(detected_count >= temp_threshold) {
  # Legge tutti i file *_operons.txt presenti nelle directory dei GSE
  all_data <- NULL
  for(gse in gse_dirs) {
    op_files <- list.files(path = gse, pattern = "*_operons.txt$", full.names = TRUE)
    if(length(op_files) > 0) {
      df <- read.table(op_files[1], header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)
      if(is.null(all_data)) {
        all_data <- df
      } else {
        all_data <- rbind(all_data, df)
      }
    }
  }

  if(!is.null(all_data)) {
    library(dplyr)
    all_data <- all_data %>%
      mutate(row_str = paste(`Start`, `Stop`, `Strand`, `Number of Genes`, Genes, sep = "\t"))

    mode_data <- all_data %>%
      group_by(Genes, row_str) %>%
      summarise(freq = n(), .groups = "drop") %>%
      group_by(Genes) %>%
      slice_max(freq, n = 1, with_ties = FALSE) %>%
      ungroup()

    mode_data <- mode_data %>%
      mutate(Start = sapply(strsplit(row_str, "\t"), `[`, 1),
             Stop = sapply(strsplit(row_str, "\t"), `[`, 2),
             Strand = sapply(strsplit(row_str, "\t"), `[`, 3),
             `Number of Genes` = sapply(strsplit(row_str, "\t"), `[`, 4)) %>%
      select(Start, Stop, Strand, `Number of Genes`, Genes)

    write.table(mode_data, file = mode_filename, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Mode operons file created:", mode_filename, "\n")
  } else {
    cat("No operons data found to generate mode file.\n")
  }
}

  if(file.exists(file.path(GSEdirectory, "_operons.txt"))) {
      new_filename <- file.path(m.o.subpath, paste0(GSE, " ", m.o.correct, "_operons.txt"))
      file.copy(file.path(GSEdirectory, "_operons.txt"), new_filename, overwrite = TRUE)
      cat("Operons file copied to:", new_filename, "\n")
    token <- TRUE
  } else {
    token <- FALSE
  }

  return(token)
}
