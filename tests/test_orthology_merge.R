args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "."
source(file.path(root, "bin/orthology_stages.R"))
base <- cbind(gene = c("a", "b"), locus = c("REF1", "REF2"), old = c("x", "y"), synonym = NA_character_)
plan <- data.frame(ordinal = 1:2, key = c("genome_000001", "genome_000002"))
one <- list(key = plan$key[1], ordinal = 1L, loci = base[, 2],
            mapping = cbind(c("A1", "A2"), c("oldA1", "oldA2"), NA_character_))
two <- list(key = plan$key[2], ordinal = 2L, loci = base[, 2],
            mapping = cbind(c("B1", "B2"), c("oldB1", "oldB2"), NA_character_))
stopifnot(identical(merge_orthology_mappings(base, plan, list(two, one)),
                    cbind(base, one$mapping, two$mapping)),
          identical(merge_orthology_mappings(base, plan[FALSE, ], list()), base))
fails <- function(pairs) inherits(tryCatch(merge_orthology_mappings(base, plan, pairs), error = identity), "error")
stopifnot(fails(list(one)), fails(list(one, one, two)))
two$loci <- rev(two$loci)
stopifnot(fails(list(one, two)))
cat("PASS: orthology merge preserves order and rejects missing/duplicate/misaligned pairs\n")
