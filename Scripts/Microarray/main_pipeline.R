#!/usr/bin/env Rscript

#############################################################################################
#################################### MAIN PIPELINE SCRIPT ###################################
#############################################################################################

# ─── 0. CATCH UNEXPECTED USER-DEFINED PARAMS ─────────────────────────────────────────────
env_names <- names(Sys.getenv())
# Recognized parameter names
allowed_env_names <- c(
  "THREADS", "MEMORY", "CORRECTION_DE", "PVALUE",
  "OPERONS_DETECTION", "OPERONS_THRESHOLD",
  "AMINO_DIFF", "SIMILARITY", "MULTITAXA", "ASMcode"
)
# Look for any env var that starts like one of our params but isn't exactly allowed
param_prefixes <- c("THREAD", "MEMORY", "CORRECTION", "PVALUE", "OPERON", "AMINO", "SIMILAR", "MULTI", "ASM")
matching_vars <- grep(paste0("^(", paste(param_prefixes, collapse = "|"), ")"), env_names, value = TRUE)
invalid_vars  <- setdiff(matching_vars, allowed_env_names)
if (length(invalid_vars) > 0) {
  stop(
    "Unexpected environment variable(s): ",
    paste(invalid_vars, collapse = ", "),
    ". Did you mean one of: ", paste(allowed_env_names, collapse = ", "), "?"
  )
}

# ─── 1. READ ENVIRONMENT VARIABLES & SET DEFAULTS ────────────────────────────────────────
raw <- function(x) Sys.getenv(x, "")
THREADS_RAW           <- raw("THREADS")
MEMORY_RAW            <- raw("MEMORY")
CORRECTION_DE_RAW     <- raw("CORRECTION_DE")
PVALUE_RAW            <- raw("PVALUE")
OPERONS_DETECTION_RAW <- raw("OPERONS_DETECTION")
AMINO_DIFF_RAW        <- raw("AMINO_DIFF")
SIMILARITY_RAW        <- raw("SIMILARITY")
MULTITAXA_RAW         <- raw("MULTITAXA")
OPERONS_THRESH_RAW    <- raw("OPERONS_THRESHOLD")
ASMcode               <- raw("ASMcode")

THREADS           <- if (nzchar(THREADS_RAW))           as.integer(THREADS_RAW)           else 8L
MEMORY            <- if (nzchar(MEMORY_RAW))            as.integer(MEMORY_RAW)            else 4000L
PVALUE            <- if (nzchar(PVALUE_RAW))            as.numeric(PVALUE_RAW)            else 0.05
OPERONS_DETECTION <- if (nzchar(OPERONS_DETECTION_RAW)) as.logical(OPERONS_DETECTION_RAW) else TRUE
n_amin            <- if (nzchar(AMINO_DIFF_RAW))        as.integer(AMINO_DIFF_RAW)        else 3L
threshold         <- if (nzchar(SIMILARITY_RAW))        as.integer(SIMILARITY_RAW)        else 30L
multispecies      <- if (nzchar(MULTITAXA_RAW))         as.logical(MULTITAXA_RAW)         else FALSE
operon_threshold  <- if (nzchar(OPERONS_THRESH_RAW))    as.integer(OPERONS_THRESH_RAW)    else 3L

# Validate CORRECTION_DE against allowed list
valid_corrections <- c("none","holm","hochberg","hommel","bonferroni","BH","BY","fdr")
if (!nzchar(CORRECTION_DE_RAW)) {
  CORRECTION_DE <- "none"
} else if (CORRECTION_DE_RAW %in% valid_corrections) {
  CORRECTION_DE <- CORRECTION_DE_RAW
} else {
  warning(
    "Invalid CORRECTION_DE '", CORRECTION_DE_RAW,
    "'; falling back to 'none'. Valid are: ",
    paste(valid_corrections, collapse = ", "), "."
  )
  CORRECTION_DE <- "none"
}

# ─── 2. VALIDATE PARAMETERS ──────────────────────────────────────────────────────────────
if (is.na(THREADS) || THREADS < 1)                   stop("THREADS must be >= 1; provided: ", THREADS)
if (is.na(MEMORY)  || MEMORY  < 1)                   stop("MEMORY must be >= 1; provided: ", MEMORY)
if (is.na(PVALUE)  || PVALUE < 0 || PVALUE > 1)      stop("PVALUE must be between 0 and 1; provided: ", PVALUE)
if (!is.logical(OPERONS_DETECTION))                  stop("OPERONS_DETECTION must be TRUE/FALSE; provided: ", OPERONS_DETECTION)
if (is.na(n_amin)   || n_amin < 0)                   stop("AMINO_DIFF must be >= 0; provided: ", n_amin)
if (is.na(threshold)|| threshold < 0)                stop("SIMILARITY must be >= 0; provided: ", threshold)
if (!is.logical(multispecies))                       stop("MULTITAXA must be TRUE/FALSE; provided: ", multispecies)
if (is.na(operon_threshold) || operon_threshold < 1) stop("OPERONS_THRESHOLD must be >= 1; provided: ", operon_threshold)
if (!nzchar(ASMcode))                                stop("ASMcode is required. Please set ASMcode.")
if (!grepl("^ASM[0-9]+v[0-9]+$", ASMcode)) {
  stop(
    "Invalid ASMcode format: '", ASMcode,
    "'. Must be 'ASM<digits>v<digits>' (e.g. ASM2046375v1)."
  )
}

# ─── 3. SET UP DIRECTORIES ─────────────────────────────────────────────────────────────
output_dir     <- "/output";        if (!dir.exists(output_dir))     dir.create(output_dir, recursive = TRUE)
tmp_dir        <- file.path(output_dir, "tmp");        if (!dir.exists(tmp_dir))        dir.create(tmp_dir)
tmp_MA_dir     <- file.path(tmp_dir, "microarray");   if (!dir.exists(tmp_MA_dir))     dir.create(tmp_MA_dir)
tmp_genome_dir <- file.path(tmp_dir, "orthologous_genome"); if (!dir.exists(tmp_genome_dir)) dir.create(tmp_genome_dir)
results_dir    <- file.path(output_dir, "results");   if (!dir.exists(results_dir))    dir.create(results_dir)
fulltab_dir    <- file.path(results_dir, "fulltab");  if (!dir.exists(fulltab_dir))    dir.create(fulltab_dir)
path_rna       <- file.path(output_dir, "RNA_Seq", "results"); if (!dir.exists(path_rna))       dir.create(path_rna, recursive = TRUE)
path_operons   <- path_rna
orthofinder    <- "/OrthoFinder/"

# ─── 4. LOAD LIBRARIES ───────────────────────────────────────────────────────────────
message("Loading required libraries...")
suppressWarnings(suppressPackageStartupMessages({
  library(limma); library(grDevices); library(base); library(Biostrings)
  library(janitor); library(S4Vectors); library(stats); library(MultNonParam)
  library(EnvStats); library(ggplot2); library(reshape2); library(GEOquery)
  library(affy); library(affyPLM); library(Biobase); library(genefilter)
  library(BiocGenerics); library(arrayQualityMetrics); library(ecolicdf)
  library(AnnotationDbi); library(tidyverse); library(dplyr)
  library(RColorBrewer); library(car); library(nortest); library(stringi)
  library(sva); library(rtracklayer); library(seqinr); library(NLP)
  library(devtools); library(ggpubr); library(utils); library(drc)
  library(nlme); library(nplr); library(schoolmath); library(data.table)
  library(biomartr); library(R.utils); library(ape); library(annotate)
  library(rentrez); library(curl); library(readxl); library(scales)
  library(googledrive); library(BiocManager); library(Rsamtools)
  library(openxlsx); library(signal); library(stringr); library(marray)
  library(swfscMisc); library(biomaRt); library(GOSemSim); library(enrichplot)
  library(org.EcK12.eg.db); library(clusterProfiler); library(downloader)
  library(oligo)
}))
message("Libraries loaded!")

# ─── 5. EXECUTION ───────────────────────────────────────────────────────────────────────
start_time <- Sys.time()
path_input  <- "/input/"

if (!file.exists(file.path(output_dir, "orthologous_genes.txt"))) {
  source("/scripts/Microarray/ortologhi_new.R")
  source("/scripts/Microarray/dowload_file_assembly.R")
  message("Starting orthologous gene search...")
  main_ortologhi(path_input, output_dir, ASMcode,
                 tmp_genome_dir, orthofinder,
                 THREADS, n_amin, threshold)
}

if (length(dir(path_input, pattern = "Microarray")) > 0) {
  message("Starting microarray analysis...")
  source("/scripts/Microarray/pipeline_MA_bash.R")
  microarray_analysis(path_input, output_dir,
                      tmp_genome_dir, tmp_MA_dir,
                      CORRECTION_DE, PVALUE)
}

if (length(dir(path_input, pattern = "RNA_Seq")) > 0) {
  message("Starting RNA-Seq analysis...")
  source("/scripts/RNA-Seq/Pipeline_bash_RNA_Seq.R", local = TRUE)
}

source("/scripts/Microarray/interexp_analysis.R")
interGSE_analysis(
  path_input,
  tmp_MA_dir,
  path_rna,
  results_dir,
  fulltab_dir,
  output_dir,
  path_operons,
  tmp_genome_dir,
  OPERONS_DETECTION,
  ASMcode,
  multispecies
)

end_time <- Sys.time()
elapsed  <- as.numeric(difftime(end_time, start_time, units = "secs"))
message(sprintf("Pipeline completed in %.2f seconds. Goodbye!", elapsed))
