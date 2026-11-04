# Isolated pairwise OrthoFinder tasks; scientific matching remains in find_ortologhi.
write_proteome <- function(genome, destination) {
  keep <- !is.na(genome[, 5])
  writeLines(as.vector(rbind(paste0(">", genome[keep, 2]), genome[keep, 5])), destination)
}

prepare_orthology <- function() {
  source("/scripts/Microarray/dowload_file_assembly.R")
  if (!file.exists("/output/assembly_summary.txt")) {
    inputs <- list.files("/input", pattern = "xlsx$", full.names = TRUE)
    assemblies <- unique(unlist(lapply(inputs, function(p) read.xlsx(p)$reference_genome)))
    reading_testo(output_dir, tmp_genome_dir, assemblies, ASMcode)
  }
  genomes <- sort(list.files(tmp_genome_dir, pattern = "\\.gb$"))
  reference <- paste0(ASMcode, ".gb")
  if (!reference %in% genomes) stop("Missing preferred reference GenBank file: ", reference)
  shared <- "/output/shared"
  dir.create(shared)
  pref <- read_gbfile(readLines(file.path(tmp_genome_dir, reference)))
  base <- pref[match(unique(pref[, 2]), pref[, 2]), 1:4, drop = FALSE]
  saveRDS(list(pref = pref, base = base, assembly = ASMcode), file.path(shared, "reference.rds"))
  write_proteome(pref, file.path(shared, paste0(ASMcode, ".fa")))
  others <- genomes[genomes != reference]
  ids <- sprintf("genome_%06d", seq_along(others))
  dir.create("/output/pairs")
  for (i in seq_along(others)) {
    folder <- file.path("/output/pairs", ids[i])
    dir.create(folder)
    assembly <- sub("\\.gb$", "", others[i])
    other <- read_gbfile(readLines(file.path(tmp_genome_dir, others[i])))
    fasta_name <- paste0(gsub(" ", "_", assembly, fixed = TRUE), ".fa")
    if (identical(fasta_name, paste0(ASMcode, ".fa"))) stop("Proteome filename collision")
    saveRDS(list(ordinal = i, key = ids[i], assembly = assembly, genome = other), file.path(folder, "genome.rds"))
    write_proteome(other, file.path(folder, fasta_name))
  }
  write.table(data.frame(ordinal = seq_along(others), key = ids), "/output/pairs.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)
}

canonical_orthogroups <- function(groups) {
  if (ncol(groups) < 5L) stop("Expected at least two species in orthogroups")
  gene_cols <- seq.int(4L, ncol(groups))
  for (j in gene_cols) {
    groups[[j]] <- vapply(as.character(groups[[j]]), function(cell) {
      if (is.na(cell) || !nzchar(cell)) return("")
      paste(sort(trimws(strsplit(cell, ",", fixed = TRUE)[[1]]), method = "radix"), collapse = ", ")
    }, character(1), USE.NAMES = FALSE)
  }
  # Stable candidate order also fixes greedy matching when alignment scores tie.
  if (nrow(groups)) groups <- groups[do.call(order, c(unname(groups[gene_cols]), list(method = "radix"))), , drop = FALSE]
  rownames(groups) <- NULL
  groups
}

run_orthology_pair <- function() {
  reference <- readRDS("/output/shared/reference.rds")
  other <- readRDS("/input/genome.rds")
  proteomes <- "/output/proteomes"
  dir.create(proteomes)
  sources <- c(list.files("/output/shared", pattern = "fa$", full.names = TRUE),
               list.files("/input", pattern = "fa$", full.names = TRUE))
  if (length(sources) != 2L || !all(file.copy(sources, proteomes))) stop("Expected two proteomes")
  status <- system2("bash", c("/adapters/run_orthofinder.sh", "-t", THREADS, "-f", shQuote(proteomes)))
  if (status != 0L) stop("OrthoFinder failed for ", other$assembly, " (exit ", status, ")")
  tables <- list.files(proteomes, pattern = "^N0\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(tables) != 1L) stop("Missing or ambiguous OrthoFinder N0.tsv")
  groups <- read.delim(tables, sep = "\t", header = TRUE, fill = TRUE)
  groups <- canonical_orthogroups(groups)
  write.table(groups, "/output/N0.canonical.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  mapped <- find_ortologhi(groups, reference$pref, other$genome, reference$assembly,
                          other$assembly, reference$base, n_amin, threshold)
  # Keep only the new columns; merge restores reference columns exactly once.
  mapping <- mapped[, seq_len(ncol(mapped)) > 4L, drop = FALSE]
  saveRDS(list(key = other$key, ordinal = other$ordinal, loci = reference$base[, 2],
               mapping = mapping), "/output/mapping.rds")
  jsonlite::write_json(list(key = other$key, ordinal = other$ordinal,
                           distance_order = "lexicographic_taxon_id_v1",
                           candidate_order = "lexicographic_gene_id_v1"), "/output/pair.json", auto_unbox = TRUE)
  file.copy(tables, "/output/N0.tsv")
  unlink(c(proteomes, "/output/shared"), recursive = TRUE)
}

merge_orthology_mappings <- function(base, plan, mappings) {
  keys <- vapply(mappings, function(x) x$key, character(1))
  if (anyDuplicated(keys) || !setequal(keys, plan$key)) stop("Missing or duplicate orthology pair results")
  result <- base
  for (key in plan$key) {
    pair <- mappings[[match(key, keys)]]
    if (pair$ordinal != plan$ordinal[match(key, plan$key)] ||
        !identical(pair$loci, base[, 2]) || nrow(pair$mapping) != nrow(result)) {
      stop("Orthology pair metadata/row mismatch: ", key)
    }
    result <- cbind(result, pair$mapping)
  }
  result
}

merge_orthology <- function() {
  reference <- readRDS("/output/shared/reference.rds")
  plan <- read.delim("/output/pairs.tsv", stringsAsFactors = FALSE)
  files <- list.files("/output/mappings", pattern = "rds$", full.names = TRUE)
  mappings <- lapply(files, readRDS)
  result <- merge_orthology_mappings(reference$base, plan, mappings)
  result <- remove_empty(result, "cols")
  write.table(result, "/output/orthologous_genes.txt", sep = "\t", quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  unlink(c("/output/shared", "/output/mappings"), recursive = TRUE)
}
