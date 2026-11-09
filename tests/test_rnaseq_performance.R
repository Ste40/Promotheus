args <- commandArgs(trailingOnly=TRUE)
root_dir <- if (length(args) > 0) args[[1]] else "."
source(file.path(root_dir, "Scripts", "RNA-Seq", "Data_retrieval.R"))
source(file.path(root_dir, "Scripts", "RNA-Seq", "BAM_parsing.R"))
source(file.path(root_dir, "Scripts", "RNA-Seq", "NDE_detection.R"))

set.seed(42)
expression_matrix <- matrix(
  sample(c(0, 0, 1:20), 6000, replace=TRUE),
  nrow=1000,
  ncol=6
)
groups <- c(1, 1, 1, 2, 2, 2)

legacy_means <- apply(expression_matrix, 1, mean)
legacy_sd <- apply(expression_matrix, 1, sd)
legacy_cv <- legacy_sd / legacy_means
legacy_cv[is.na(legacy_cv)] <- 0
statistics <- expression_statistics(expression_matrix)
stopifnot(isTRUE(all.equal(statistics$means, legacy_means, tolerance=1e-12)))
stopifnot(isTRUE(all.equal(statistics$sd, legacy_sd, tolerance=0)))
stopifnot(isTRUE(all.equal(statistics$cv, legacy_cv, tolerance=1e-12)))

legacy_kruskal <- apply(expression_matrix, 1, function(values) {
  kruskal.test(values ~ groups)$p.value
})
optimized_kruskal <- kruskal_pvalues(expression_matrix, groups)
stopifnot(isTRUE(all.equal(optimized_kruskal, legacy_kruskal, tolerance=1e-15)))

group_length_error <- tryCatch(
  kruskal_pvalues(expression_matrix, groups[-1]),
  error = identity
)
single_group_error <- tryCatch(
  kruskal_pvalues(expression_matrix, rep(1, ncol(expression_matrix))),
  error = identity
)
stopifnot(
  inherits(group_length_error, "error"),
  inherits(single_group_error, "error")
)

legacy_percentiles <- apply(expression_matrix, 2, function(values) {
  distribution <- ecdf(sort(values))
  vapply(values, distribution, numeric(1))
})
optimized_percentiles <- expression_percentiles(expression_matrix)
stopifnot(isTRUE(all.equal(optimized_percentiles, legacy_percentiles, tolerance=0)))

samtools_arguments <- list()
fake_samtools <- function(command, args) {
  stopifnot(command == "samtools")
  samtools_arguments[[length(samtools_arguments) + 1L]] <<- args
  0L
}
stopifnot(run_samtools(c("sort", "-@", "8", "input.bam"), runner=fake_samtools))
stopifnot("-@" %in% samtools_arguments[[1]])
stopifnot("8" %in% samtools_arguments[[1]])

infer_arguments <- NULL
fake_infer <- function(command, args, stdout, stderr) {
  stopifnot(command == "infer_experiment.py", isTRUE(stdout), isTRUE(stderr))
  infer_arguments <<- args
  "This is PairEnd Data"
}
infer_output <- infer_experiment_stats("sample.sorted.bam", runner=fake_infer)
stopifnot(
  identical(infer_output, "This is PairEnd Data"),
  identical(
    infer_arguments,
    c("-s", "5000000", "-i", shQuote("sample.sorted.bam"), "-r", shQuote("gff.bed"))
  )
)

gzip_arguments <- list()
fake_gzip <- function(command, args) {
  stopifnot(command == "gzip")
  gzip_arguments[[length(gzip_arguments) + 1L]] <<- args
  0L
}
stopifnot(compress_fastq_files(c("one.fastq", "two.fastq"), 1, runner=fake_gzip))
stopifnot(length(gzip_arguments) == 2L)
stopifnot(all(vapply(
  gzip_arguments,
  function(arguments) all(c("-1", "-f") %in% arguments),
  logical(1)
)))

compression_dir <- tempfile("promotheus-compression-")
dir.create(compression_dir)
on.exit(unlink(compression_dir, recursive=TRUE, force=TRUE), add=TRUE)
fastq_paths <- file.path(compression_dir, c("pair_1.fastq", "pair_2.fastq"))
for (path in fastq_paths) {
  writeLines(c("@read", "ACGT", "+", "!!!!"), path)
}
stopifnot(compress_fastq_files(fastq_paths, 2))
stopifnot(all(file.exists(paste0(fastq_paths, ".gz"))))
stopifnot(all(file.info(paste0(fastq_paths, ".gz"))$size > 0))

cat("PASS: optimized RNA-Seq calculations preserve legacy results\n")
