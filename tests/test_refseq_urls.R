args <- commandArgs(trailingOnly = TRUE)
root_dir <- args[[1]]

scripts <- c(
  file.path(root_dir, "Scripts", "Microarray", "dowload_file_assembly.R"),
  file.path(root_dir, "Scripts", "RNA-Seq", "Data_retrieval.R"),
  file.path(root_dir, "Scripts", "Microarray", "interexp_analysis.R")
)

ftp_path <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/176/835/",
  "GCF_003176835.1_ASM317683v1"
)
expected <- paste0(
  ftp_path,
  "/GCF_003176835.1_ASM317683v1_genomic.gbff.gz"
)

for (script in scripts) {
  environment <- new.env(parent = globalenv())
  sys.source(script, envir = environment)

  without_slash <- environment$refseq_file_location(ftp_path, "_genomic.gbff.gz")
  with_slash <- environment$refseq_file_location(paste0(ftp_path, "/"), "_genomic.gbff.gz")

  stopifnot(
    identical(without_slash$url, expected),
    identical(with_slash$url, expected),
    identical(with_slash$filename, "GCF_003176835.1_ASM317683v1_genomic.gbff.gz")
  )

  invalid_path_fails <- inherits(
    try(environment$refseq_file_location("na", "_genomic.gbff.gz"), silent = TRUE),
    "try-error"
  )
  stopifnot(invalid_path_fails)
}

message("PASS: RefSeq URL compatibility regression test")
