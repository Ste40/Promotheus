args <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(args) > 0L) args[[1]] else "."

source(file.path(root_dir, "Scripts", "RNA-Seq", "NDE_detection.R"))
source(file.path(root_dir, "Scripts", "Microarray", "interexp_analysis.R"))

rna_result <- data.frame(
  GeneID = c("gene_a", "gene_b"),
  locus_tag = c("LOC_A", "LOC_B"),
  sample_1 = c(1, 2),
  means = c(1, 2),
  sd = c(0.1, 0.2),
  diff_signif = c("NO", "YES"),
  stringsAsFactors = FALSE
)

annotated <- add_operon_status_column(rna_result, "OPERON_ANALYSIS_DISABLED")
stopifnot(
  nrow(annotated) == nrow(rna_result),
  all(colnames(rna_result) %in% colnames(annotated)),
  identical(annotated$GeneID, rna_result$GeneID),
  identical(annotated$locus_tag, rna_result$locus_tag),
  all(annotated$n.genes == "OPERON_ANALYSIS_DISABLED")
)

temporary_directory <- tempfile("promotheus-operons-disabled-")
dir.create(temporary_directory)
on.exit(unlink(temporary_directory, recursive = TRUE), add = TRUE)

valid_path <- file.path(temporary_directory, "valid_full_tab.txt")
write.table(rna_result, valid_path, sep = "\t", row.names = FALSE, quote = FALSE)
validated <- read_expression_full_tab(valid_path)
stopifnot(
  identical(validated$locus_tag, rna_result$locus_tag),
  is_valid_rnaseq_full_tab(valid_path)
)

invalid_path <- file.path(temporary_directory, "invalid_full_tab.txt")
write.table(
  data.frame(n.genes = rep("WAITING_FOR_OPERONS_TABLE", 2)),
  invalid_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

validation_error <- tryCatch(
  {
    read_expression_full_tab(invalid_path)
    NULL
  },
  error = identity
)
stopifnot(
  inherits(validation_error, "error"),
  grepl("missing columns: locus_tag, means, diff_signif", conditionMessage(validation_error), fixed = TRUE),
  !is_valid_rnaseq_full_tab(invalid_path),
  !is_valid_rnaseq_full_tab(file.path(temporary_directory, "missing_full_tab.txt"))
)

cat("PASS: operons-disabled RNA-Seq full_tab regression tests\n")
