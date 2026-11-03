args <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(args) > 0L) args[[1]] else "."
source(file.path(root_dir, "Scripts", "Microarray", "interexp_analysis.R"))

expression_matrix <- matrix(
  c(2, 4, 8, 3, 6, 9, 5, 5, 10),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(NULL, c("GSE1", "GSE2", "GSE3"))
)
complete_table <- cbind(
  Geneid = c("gene_a", "gene_b", "gene_c"),
  locus_comuni = c("LOC_A", "LOC_B", "LOC_C"),
  expression_matrix,
  DE = c("NO", "YES", "NO")
)

legacy_means <- apply(expression_matrix, 1, mean)
legacy_sd <- apply(expression_matrix, 1, sd)
legacy_percentiles <- apply(expression_matrix, 2, function(values) {
  rank(values) / nrow(expression_matrix) * 100
})
legacy_percentile_means <- apply(legacy_percentiles, 1, mean)
legacy_percentile_sd <- apply(legacy_percentiles, 1, sd)

result <- create_final_table(expression_matrix, complete_table)
stopifnot(
  identical(result$Gene_name, complete_table[, 1]),
  identical(result$GeneID, complete_table[, 2]),
  isTRUE(all.equal(as.numeric(result$Mean_Expression), legacy_means, tolerance = 1e-12)),
  isTRUE(all.equal(as.numeric(result$CV), legacy_sd, tolerance = 1e-12)),
  isTRUE(all.equal(as.numeric(result$Mean_Percentile), legacy_percentile_means, tolerance = 1e-12)),
  isTRUE(all.equal(as.numeric(result$Percentile_SD), legacy_percentile_sd, tolerance = 1e-12)),
  identical(result$DE, complete_table[, ncol(complete_table)])
)

combat_calls <- 0L
fake_combat <- function(dat, batch, mod, par.prior, prior.plots, mean.only) {
  combat_calls <<- combat_calls + 1L
  stopifnot(
    is.null(mod),
    isTRUE(par.prior),
    identical(prior.plots, FALSE),
    isTRUE(mean.only)
  )
  dat + 1
}

single_batch <- adjust_batch_effect(expression_matrix, rep("GSE1", ncol(expression_matrix)), fake_combat)
multi_batch <- adjust_batch_effect(expression_matrix, c("GSE1", "GSE1", "GSE2"), fake_combat)
length_error <- tryCatch(
  adjust_batch_effect(expression_matrix, c("GSE1", "GSE2"), fake_combat),
  error = identity
)
stopifnot(
  identical(single_batch, expression_matrix),
  identical(multi_batch, expression_matrix + 1),
  combat_calls == 1L,
  inherits(length_error, "error")
)

fasta <- strsplit("ACGTACGTAA", "", fixed = TRUE)[[1]]
stopifnot(
  identical(promoter_sequence(fasta, 3L, 5L, "+", 4L), "AAAC"),
  identical(promoter_sequence(fasta, 3L, 9L, "-", 4L), "AACG"),
  nchar(promoter_sequence(fasta, 1L, 1L, "+")) == 300L
)

# ape::read.gff uses named columns, including attributes rather than V9.
# A prefix match must not confuse LOC_1 with LOC_10 or an old locus tag.
gff <- data.frame(attributes = c("ID=gene-LOC_10;locus_tag=LOC_10;old_locus_tag=LOC_1",
                               "ID=gene-LOC_1;locus_tag=LOC_1",
                               "ID=cds-LOC_1;Parent=gene-LOC_1;locus_tag=LOC_1"))
stopifnot(identical(reference_gene_rows(gff, c("LOC_1", "LOC_10")), c(2L, 1L)))
names(gff) <- "V9"
stopifnot(identical(reference_gene_rows(gff, "LOC_1"), 2L),
          inherits(tryCatch(reference_gene_rows(gff, "MISSING"), error = identity), "error"))

cat("PASS: inter-experiment calculations and GFF reference matching\n")
