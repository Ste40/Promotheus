args <- commandArgs(trailingOnly=TRUE)
root_dir <- if (length(args) > 0) args[[1]] else "."
source(file.path(root_dir, "Scripts", "RNA-Seq", "Data_retrieval.R"))

test_dir <- tempfile("promotheus-fastq-retry-")
dir.create(test_dir, recursive=TRUE)
on.exit(unlink(test_dir, recursive=TRUE, force=TRUE), add=TRUE)

attempts <- 0L
successful_runner <- function(command, args) {
  stopifnot(command == "fasterq-dump")
  attempts <<- attempts + 1L
  if (attempts == 3L) {
    writeLines("@read\nACGT\n+\n!!!!", file.path(test_dir, "SRR000001.fastq"))
    return(0L)
  }
  writeLines("partial", file.path(test_dir, "SRR000001.fastq"))
  3L
}

stopifnot(download_fastq_with_retry(
  "SRR000001", test_dir, 2,
  max_attempts=3L, base_delay=0,
  runner=successful_runner, sleep_fn=function(seconds) NULL
))
stopifnot(attempts == 3L)
stopifnot(length(fastq_files_for_accession("SRR000001", test_dir)) == 1L)

unlink(fastq_files_for_accession("SRR000001", test_dir), force=TRUE)
attempts <- 0L
failing_runner <- function(command, args) {
  attempts <<- attempts + 1L
  writeLines("partial", file.path(test_dir, "SRR000002.fastq"))
  3L
}

stopifnot(!download_fastq_with_retry(
  "SRR000002", test_dir, 2,
  max_attempts=3L, base_delay=0,
  runner=failing_runner, sleep_fn=function(seconds) NULL
))
stopifnot(attempts == 3L)
stopifnot(length(fastq_files_for_accession("SRR000002", test_dir)) == 0L)

writeLines("compressed", file.path(test_dir, "SRR000003.fastq.gz"))
writeLines("compressed", file.path(test_dir, "SRR000004_1.fastq.gz"))
writeLines("compressed", file.path(test_dir, "SRR000004_2.fastq.gz"))
stopifnot(length(fastq_files_for_accession("SRR000003", test_dir, compressed=TRUE)) == 1L)
stopifnot(length(fastq_files_for_accession("SRR000004", test_dir, compressed=TRUE)) == 2L)
stopifnot(fastq_accession_complete("SRR000003", test_dir, compressed=TRUE))
stopifnot(fastq_accession_complete("SRR000004", test_dir, compressed=TRUE))

invisible(file.create(file.path(test_dir, "SRR000005.fastq.gz")))
stopifnot(!fastq_accession_complete("SRR000005", test_dir, compressed=TRUE))

for (name in c("SRR000006.fastq", "SRR000006_1.fastq", "SRR000006_2.fastq")) {
  writeLines("reads", file.path(test_dir, name))
}
normalized_pair <- normalize_split_fastq_files("SRR000006", test_dir)
stopifnot(
  identical(sort(basename(normalized_pair)), c("SRR000006_1.fastq", "SRR000006_2.fastq")),
  !file.exists(file.path(test_dir, "SRR000006.fastq")),
  fastq_accession_complete("SRR000006", test_dir, compressed=FALSE)
)

cat("PASS: fasterq-dump retry and partial-file cleanup tests\n")
