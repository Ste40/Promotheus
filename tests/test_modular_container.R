# Read-only environment/adapter contract smoke test; no NCBI requests.
Sys.setenv(ASMcode = "ASM317683v1", THREADS = "1", MEMORY = "4000")
main_lines <- readLines("/scripts/Microarray/main_pipeline.R")
boundary <- grep("^start_time <-", main_lines)
stopifnot(length(boundary) == 1L)
eval(parse(text = main_lines[seq_len(boundary - 1L)]), envir = globalenv())
suppressPackageStartupMessages({ library(Rsubread); library(DESeq2); library(ProGenome) })
for (script in c("Data_retrieval.R", "MetaData_retrieval.R", "BAM_parsing.R",
                 "Operons_detection.R", "NDE_detection.R", "merge_operons_table.R")) {
  source(file.path("/scripts/RNA-Seq", script))
}
source("/scripts/Microarray/pipeline_MA_bash.R")
source("/scripts/Microarray/interexp_analysis.R")
source("/scripts/Microarray/ortologhi_new.R")
invisible(parse("/adapters/run_stage.R"))
input <- read.xlsx("/input/RNA_Seq.xlsx")
stopifnot(nrow(input) == 1L,
          all(c("Accession", "GSE", "taxon", "reference_genome") %in% names(input)),
          length(input$Accession) == 1L,
          identical(as.character(input$Accession), "GSE000002"),
          identical(paste0("GSE", input$GSE), as.character(input$Accession)))
for (name in c("main_ortologhi", "microarray_analysis", "MetaData_retrieval",
               "Data_retrieval", "BAM_parsing", "Operons_detection", "NDE_detection",
               "merge_operons", "interGSE_analysis")) {
  stopifnot(is.function(get(name)))
}
cat("PASS: modular adapters, rewritten workbook and scientific packages load in pinned image\n")
