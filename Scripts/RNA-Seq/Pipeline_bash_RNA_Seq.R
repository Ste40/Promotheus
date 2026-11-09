

#######################################
# Welcome to Ste40cm RNA-Seq Pipeline #
#######################################


INPUT_XLSX <- list.files(path = "/input", pattern = "RNA_Seq", full.names = TRUE)

setwd("/scripts/RNA-Seq") # Pipeline directory

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

if (!isTRUE(Packages_installation(packages))) {
  stop("Oooops! Something went wrong during package installation")
} else {
  message("Installation successful!")
}

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

microorganism_table <- read.xlsx(INPUT_XLSX)

xlsx <- data.frame(m.o.name = character(), GSE = character(), Experimental_condition = character(), Groups = character(),
  Accession_codes = character())


if (!file.exists(file.path("/output/RNA_Seq/Experimental_Groups.xlsx"))) {
  message("I'm retrieving metadata, please wait")
  for (g in seq_len(nrow(microorganism_table))) {
    xlsx_app <- suppressWarnings(MetaData_retrieval(microorganism_table, g, xlsx))
    xlsx <- rbind(xlsx, xlsx_app)
  }

  write.xlsx(xlsx, file.path("/output/RNA_Seq/Experimental_Groups.xlsx"), overwrite = TRUE)
  if (nrow(xlsx) != length(str_split(microorganism_table$GSE, ";")[])) {
    message("Oooops! Something went wrong during metadata retrieval, please check 'Experimental_Groups.xlsx' file")
  } else {
    message("Success: metadata retrieved!")
  }
} else {
  message("Success: Experimental Groups already defined!")
}


cat(
  "Before proceeding, check if the groups and replicates have been assigned correctly ",
  "within table \"Experimental_Groups.xlsx\".\n",
  "If not, please edit the file and overwrite it. When you are done, press enter:"
)
invisible(scan("stdin", character(), nlines = 1, quiet = TRUE))


#######################################


Conditions <- read.xlsx(list.files(path = "/output/RNA_Seq", pattern = "*Groups.xlsx", full.names = TRUE)[1])
GSEs <- Conditions[, 2]
gse_count <- aggregate(GSE ~ m.o.name, data = Conditions, FUN = function(x) length(unique(x)))
original_threshold <- operon_threshold
od_full_tabs <- file.path("/output/RNA_Seq/results")
m.o.correct.all <- str_replace_all(Conditions$m.o.name, "[[:punct:]]", " ")
full_tab_paths <- file.path(od_full_tabs, paste(GSEs, m.o.correct.all, "full_tab.txt"))

max_gse_parallel <- suppressWarnings(as.integer(Sys.getenv("MAX_GSE_PARALLEL", "1")))
if (is.na(max_gse_parallel) || max_gse_parallel < 1L) {
  stop("MAX_GSE_PARALLEL must be a positive integer")
}

# GSEs of the same organism remain serial because operon-mode aggregation depends on their order.
# Different organisms use separate directories and can therefore execute safely in parallel.
organism_groups <- split(seq_along(GSEs), m.o.correct.all)
parallel_workers <- min(max_gse_parallel, length(organism_groups), as.integer(THREADS))
threads_per_worker <- max(1L, floor(as.integer(THREADS) / parallel_workers))
memory_per_worker <- max(512L, floor(as.integer(MEMORY) / parallel_workers))

process_rnaseq_gse <- function(r, threads, memory) {
  temp_threshold <- original_threshold
  GSE <- GSEs[r]
  Name_ID <- Conditions$m.o.name[r]
  Assembly_v <- microorganism_table$reference_genome[r]
  m.o.correct <- m.o.correct.all[r]
  GSEconditions <- Conditions[r, c("Experimental_condition", "Groups")]
  m.o.subpath <- file.path("/output/RNA_Seq/", m.o.correct)
  Acclist <- str_split(Conditions$Accession_codes[r], ",")[[1]]
  GSEdirectory <- file.path(m.o.subpath, paste(GSE, m.o.correct))
  full_tab_path <- full_tab_paths[r]

  if (is_valid_rnaseq_full_tab(full_tab_path)) {
    message(paste(GSE, "has already been analyzed; proceeding with the next one.."))
    return(TRUE)
  }
  if (file.exists(full_tab_path)) {
    warning("Existing RNA-Seq result for ", GSE, " is incomplete or malformed and will be regenerated: ",
      full_tab_path)
  }

  suppressWarnings(dir.create(m.o.subpath, recursive = TRUE))
  message(paste("I'm downloading the required files for", Name_ID, GSE))
  retr <- Data_retrieval(GSE, Acclist, GSEdirectory, m.o.subpath, args, threads, GSEconditions,
    Assembly_v, Name_ID)
  if (identical(retr, "skip")) {
    return(FALSE)
  }
  if (isTRUE(retr)) {
    message("Success: files downloaded!")
  }

  message(paste("I'm aligning the reads for", GSE, "to reference genome, please wait"))
  strandness <- BAM_parsing(GSE, m.o.subpath, GSEdirectory, Acclist, threads)
  if (strandness[[1]][1] != "This is SingleEnd Data" && strandness[[1]][1] != "This is PairEnd Data") {
    stop("Oooops! Something went wrong during alignment or .bam parsing")
  }
  message("Success: .BAM file analyzed!")

  if (isTRUE(OPERONS_DETECTION)) {
    message("I'm detecting operons, please wait")
    if (!isTRUE(Operons_detection(m.o.subpath, threads, memory, Acclist, GSEdirectory, GSEconditions,
      temp_threshold, gse_count, GSE, m.o.correct))) {
      warning(paste("Oooops! Something went wrong during Operon detection for", basename(GSEdirectory),
        "skipping operons detection for this GSE"))
    } else {
      message("Success: Operons detected!")
    }
  } else {
    message("OPERONS_DETECTION is set to FALSE, operons detection will be skipped")
  }

  message("I'm detecting Not-Differentially Expressed (NDE) genes, please wait")
  if (!isTRUE(NDE_detection(GSE, GSEdirectory, strandness, threads, m.o.subpath, GSEconditions, m.o.correct,
    od_full_tabs, CORRECTION_DE, OPERONS_DETECTION, PVALUE, Acclist))) {
    stop("Oooops! Something went wrong during NDE identification")
  }
  message("Success: NDE genes detected!")

  if (isTRUE(OPERONS_DETECTION)) {
    tryCatch({
      message("I'm updating the GSE table based on operons information, please wait")
      merge_operons(m.o.subpath, od_full_tabs, GSEdirectory, GSE, m.o.correct)
      message("Success: operon annotations updated!")
    }, error = function(e) {
      stop("Oooops! Something went wrong during operons merge: ", e$message)
    })
  }

  message(paste(GSE, "analysis completed!"))
  TRUE
}

process_organism_group <- function(indices) {
  vapply(indices, process_rnaseq_gse, logical(1), threads = threads_per_worker, memory = memory_per_worker)
}

if (parallel_workers > 1L && .Platform$OS.type != "windows") {
  message("Processing ", length(GSEs), " RNA-Seq GSEs in ", parallel_workers,
    " organism-isolated workers (", threads_per_worker, " threads each)")
  group_results <- parallel::mclapply(
    organism_groups,
    process_organism_group,
    mc.cores = parallel_workers,
    mc.preschedule = FALSE
  )
} else {
  group_results <- lapply(organism_groups, process_organism_group)
}

worker_errors <- vapply(group_results, inherits, logical(1), "try-error")
if (any(worker_errors)) {
  stop("One or more RNA-Seq organism workers failed: ", paste(names(group_results)[worker_errors], collapse = ", "))
}

valid_full_tabs <- vapply(full_tab_paths, is_valid_rnaseq_full_tab, logical(1))
if (all(valid_full_tabs)) {
  message("All the GSEs have been analyzed!")
} else {
  stop(
    "RNA-Seq datasets without a valid full_tab result: ",
    paste(GSEs[!valid_full_tabs], collapse = ", ")
  )
}
