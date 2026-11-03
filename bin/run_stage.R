#!/usr/bin/env Rscript
# Task-local adapters around the existing scientific functions. No shared output mount.
stage <- commandArgs(trailingOnly = TRUE)[[1]]
main_lines <- readLines("/scripts/Microarray/main_pipeline.R")
boundary <- grep("^start_time <-", main_lines)
stopifnot(length(boundary) == 1L)
eval(parse(text = main_lines[seq_len(boundary - 1L)]), envir = globalenv())
path_input <- "/input/"
setwd("/output")

if (stage %in% c("orthology_prepare", "orthology_pair", "orthology_merge")) {
  source("/scripts/Microarray/ortologhi_new.R")
  source("/adapters/orthology_stages.R")
  switch(stage, orthology_prepare = prepare_orthology(), orthology_pair = run_orthology_pair(),
         orthology_merge = merge_orthology())
} else if (stage == "orthology") {
  source("/scripts/Microarray/ortologhi_new.R")
  main_ortologhi(path_input, output_dir, ASMcode, tmp_genome_dir, orthofinder,
                THREADS, n_amin, threshold)
} else if (stage == "microarray") {
  source("/scripts/Microarray/pipeline_MA_bash.R")
  microarray_analysis(path_input, output_dir, tmp_genome_dir, tmp_MA_dir, CORRECTION_DE, PVALUE)
  tabs <- list.files(tmp_MA_dir, pattern = "_full_tab.txt$", full.names = TRUE)
  if (length(tabs) != 1L) stop("Microarray task did not produce its requested GSE full_tab")
  file.copy(tabs, fulltab_dir, overwrite = TRUE)
} else if (stage %in% c("rna_download", "rna_align", "rna_analyze")) {
  suppressPackageStartupMessages({
    library(Rsubread); library(DESeq2); library(ProGenome)
  })
  for (script in c("Data_retrieval.R", "MetaData_retrieval.R", "BAM_parsing.R",
                   "Operons_detection.R", "NDE_detection.R", "merge_operons_table.R")) {
    source(file.path("/scripts/RNA-Seq", script))
  }
  if (stage == "rna_download") {
    input <- read.xlsx(list.files(path_input, pattern = "RNA_Seq.*xlsx$", full.names = TRUE))
    stopifnot(nrow(input) == 1L,
              all(c("Accession", "GSE", "taxon", "reference_genome") %in% names(input)))
    conditions <- MetaData_retrieval(input, 1L, data.frame())
    stopifnot(nrow(conditions) == 1L, length(conditions$GSE) == 1L,
              identical(as.character(conditions$GSE), as.character(input$Accession)))
    organism <- str_replace_all(conditions$m.o.name, "[[:punct:]]", " ")
    gse <- conditions$GSE[[1]]
    organism_path <- file.path("/output/RNA_Seq", organism)
    gse_path <- file.path(organism_path, paste(gse, organism))
    dir.create(organism_path, recursive = TRUE, showWarnings = FALSE)
    accessions <- str_split(conditions$Accession_codes, ",")[[1]]
    groups <- conditions[, c("Experimental_condition", "Groups"), drop = FALSE]
    if (!isTRUE(Data_retrieval(gse, accessions, gse_path, organism_path, character(),
                              THREADS, groups, input$reference_genome[[1]], input$taxon[[1]]))) {
      stop("RNA-Seq download incomplete for ", gse)
    }
    saveRDS(list(conditions = conditions, input = input, organism = organism,
                 gse = gse, accessions = accessions), "/output/state.rds")
    write.xlsx(conditions, "/output/RNA_Seq/Experimental_Groups.xlsx", overwrite = TRUE)
  } else if (stage == "rna_align") {
    state <- readRDS("/output/state.rds")
    organism_path <- file.path("/output/RNA_Seq", state$organism)
    gse_path <- file.path(organism_path, paste(state$gse, state$organism))
    state$strandness <- BAM_parsing(state$gse, organism_path, gse_path, state$accessions, THREADS)
    if (!(state$strandness[[1]][1] %in% c("This is SingleEnd Data", "This is PairEnd Data"))) {
      stop("Alignment/strand inference failed for ", state$gse)
    }
    saveRDS(state, "/output/state.rds")
  } else {
    # Python stages states in original input order. Operon consensus retains this order.
    state_files <- list.files("/output/states", pattern = "rds$", full.names = TRUE)
    states <- lapply(state_files, readRDS)
    conditions <- do.call(rbind, lapply(states, function(x) x$conditions))
    gse_count <- aggregate(GSE ~ m.o.name, data = conditions, FUN = function(x) length(unique(x)))
    write.xlsx(conditions, "/output/RNA_Seq/Experimental_Groups.xlsx", overwrite = TRUE)
    for (state in states) {
      organism_path <- file.path("/output/RNA_Seq", state$organism)
      gse_path <- file.path(organism_path, paste(state$gse, state$organism))
      groups <- state$conditions[, c("Experimental_condition", "Groups"), drop = FALSE]
      if (isTRUE(OPERONS_DETECTION)) {
        detected <- Operons_detection(organism_path, THREADS, MEMORY, state$accessions,
          gse_path, groups, operon_threshold, gse_count, state$gse, state$organism)
        if (!isTRUE(detected)) warning("Operon detection incomplete for ", state$gse)
      }
      if (!isTRUE(NDE_detection(state$gse, gse_path, state$strandness, THREADS, organism_path,
          groups, state$organism, path_rna, CORRECTION_DE, OPERONS_DETECTION, PVALUE))) {
        stop("NDE analysis failed for ", state$gse)
      }
      if (isTRUE(OPERONS_DETECTION)) {
        merge_operons(organism_path, path_rna, gse_path, state$gse, state$organism)
      }
      expected <- file.path(path_rna, paste(state$gse, state$organism, "full_tab.txt"))
      if (!is_valid_rnaseq_full_tab(expected)) stop("Invalid RNA-Seq full_tab: ", expected)
    }
  }
} else if (stage == "integrate") {
  source("/scripts/Microarray/interexp_analysis.R")
  interGSE_analysis(path_input, tmp_MA_dir, path_rna, results_dir, fulltab_dir,
                   output_dir, path_rna, tmp_genome_dir, OPERONS_DETECTION, ASMcode, multispecies)
  if (!file.exists(file.path(results_dir, "final_table.xlsx"))) stop("Missing integrated workbook")
} else {
  stop("Unknown stage: ", stage)
}
