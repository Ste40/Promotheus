args <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(args) > 0L) args[[1]] else "."

run_tests <- function() {
  fixture <- tempfile("promotheus-integration-")
  dir.create(fixture)
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
  ma_dir <- file.path(fixture, "microarray")
  rna_dir <- file.path(fixture, "rnaseq")
  dir.create(ma_dir)
  dir.create(rna_dir)
  orthologues <- data.frame(gene = c("a", "b"), reference = c("REF_A", "REF_B"),
                           other = c("RNA_A", "RNA_B"))

  # Redirect only the fixed container orthology path; all expression tables are
  # written and parsed using the real producers' column layouts.
  scope <- new.env(parent = globalenv())
  scope$read.table <- function(file, ...) {
    if (identical(file, "/output/orthologous_genes.txt")) return(orthologues)
    if (startsWith(file, "/output/RNA_Seq/results/")) file <- file.path(rna_dir, basename(file))
    utils::read.table(file, ...)
  }
  scope$dir <- function(path, ...) {
    if (identical(path, "/output/RNA_Seq")) return(c("results", "Toy"))
    if (identical(path, "/output/RNA_Seq/results")) path <- rna_dir
    base::dir(path, ...)
  }
  sys.source(file.path(root_dir, "Scripts/Microarray/interexp_analysis.R"), scope)
  ma <- data.frame(Geneid = c("a", "b"), locus_tag = c("REF_A", "REF_B"),
                   GSM1 = c(10, 20), GSM2 = c(30, 40), means = c(20, 30),
                   sd = 1, diff_signif = c("NO", "YES"))
  # NDE_detection appends locus_tag after sample/statistic columns.
  rna <- data.frame(GeneID = c("a", "b"), SRR1 = c(1, 2), SRR2 = c(3, 4),
                    SRR3 = c(5, 6), SRR4 = c(7, 8), means = c(4, 5), sd = 1,
                    diff_signif = c("YES", "NO"), locus_tag = c("REF_A", "REF_B"))
  write.table(ma, file.path(ma_dir, "GSE1_full_tab.txt"), sep = "\t", row.names = FALSE)
  rna_path <- file.path(rna_dir, "GSE2 Toy full_tab.txt")
  check <- function(table, expected) {
    write.table(table, rna_path, sep = "\t", row.names = FALSE)
    result <- scope$unique_matrix(ma_dir, rna_dir, fixture, FALSE)
    aligned <- result[[1]][match(c("REF_A", "REF_B"), result[[1]][, "locus_comuni"]), , drop = FALSE]
    stopifnot(length(result[[2]]) == 6L,
              identical(unname(as.numeric(table(result[[2]]))), c(2, 4)),
              identical(unname(matrix(as.numeric(aligned[, 3:8]), nrow = 2)),
                        unname(as.matrix(expected))),
              identical(unname(aligned[, "DE"]), c("YES", "YES")))
  }
  expected <- cbind(ma[, c("GSM1", "GSM2")], rna[, c("SRR1", "SRR2", "SRR3", "SRR4")])
  check(rna, expected)
  # Reordered input rows and a different assembly must map to reference loci
  # even when every gene has an orthologue (no unmapped rows to discard).
  rna$locus_tag <- c("RNA_A", "RNA_B")
  check(rna[2:1, ], expected)
  # Also retain the two-leading-identifiers layout used by older full tables.
  check(rna[, c("GeneID", "locus_tag", "SRR1", "SRR2", "SRR3", "SRR4",
                "means", "sd", "diff_signif")], expected)
  rna$n.genes <- c("first gene in operon", "not first gene in operon")
  write.table(rna, rna_path, sep = "\t", row.names = FALSE)
  filtered <- scope$unique_matrix(ma_dir, rna_dir, fixture, TRUE)[[1]]
  stopifnot(nrow(filtered) == 1L, filtered[1, "locus_comuni"] == "REF_A",
            identical(as.numeric(filtered[1, 3:8]), as.numeric(expected[1, ])))
}
run_tests()
cat("PASS: integration retains every MA/RNA sample and maps reference loci\n")
