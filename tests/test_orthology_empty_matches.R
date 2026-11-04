args <- commandArgs(trailingOnly = TRUE)
root_dir <- args[[1]]

environment <- new.env(parent = globalenv())
environment$is_empty <- function(value) length(value) == 0
sys.source(
  file.path(root_dir, "Scripts", "Microarray", "ortologhi_new.R"),
  envir = environment
)

orthogroups <- data.frame(
  HOG = "N0.HOG0000000",
  OG = "OG0000000",
  parent = "n0",
  "ASM317683v1" = "P1, P2",
  "ASM3872123v1" = "O1",
  check.names = FALSE
)

preferred <- matrix(
  c(
    "gene1", "P1", NA, NA, "AAAA",
    "gene2", "P2", NA, NA, "CCCC"
  ),
  nrow = 2,
  byrow = TRUE
)
other <- matrix(c("gene3", "O1", NA, NA, "A"), nrow = 1)
final_table <- matrix(c("gene0", "P0", NA, NA), nrow = 1)

result <- environment$find_ortologhi(
  orthogroups,
  preferred,
  other,
  "ASM317683v1",
  "ASM3872123v1",
  final_table,
  n_amin = 0,
  threashold = 100
)

stopifnot(
  identical(dim(result), dim(final_table)),
  identical(as.character(result), as.character(final_table))
)

if (!methods::isClass("MockAlignment")) {
  methods::setClass(
    "MockAlignment",
    slots = c(score = "numeric", aligned_pattern = "character")
  )
}
environment$pairwiseAlignment <- function(prot1, prot2, type) {
  methods::new(
    "MockAlignment",
    score = 10,
    aligned_pattern = as.character(prot1)
  )
}
environment$nmatch <- function(alignment) nchar(alignment@aligned_pattern)
environment$pattern <- function(alignment) alignment@aligned_pattern

tie_orthogroups <- data.frame(
  HOG = "N0.HOG0000001",
  OG = "OG0000001",
  parent = "n0",
  "ASM317683v1" = "P1, P2, P3, P4",
  "ASM207885v1" = "O1, O2, O3",
  check.names = FALSE
)
tie_preferred <- matrix(
  c(
    "gene1", "P1", NA, NA, "AAAA",
    "gene2", "P2", NA, NA, "AAAA",
    "gene3", "P3", NA, NA, "AAAA",
    "gene4", "P4", NA, NA, "AAAA"
  ),
  nrow = 4,
  byrow = TRUE
)
tie_other <- matrix(
  c(
    "other1", "O1", NA, NA, "AAAA",
    "other2", "O2", NA, NA, "AAAA",
    "other3", "O3", NA, NA, "AAAA"
  ),
  nrow = 3,
  byrow = TRUE
)
tie_final <- tie_preferred[, 1:4, drop = FALSE]

tie_result <- environment$find_ortologhi(
  tie_orthogroups,
  tie_preferred,
  tie_other,
  "ASM317683v1",
  "ASM207885v1",
  tie_final,
  n_amin = 3,
  threashold = 30
)

stopifnot(
  nrow(tie_result) == 4,
  ncol(tie_result) == 7,
  sum(!is.na(tie_result[, 5])) == 3,
  all(na.omit(tie_result[, 5]) %in% c("O1", "O2", "O3"))
)

# Fixed OrthoFinder groups must give the same map when each genome is matched
# independently and results arrive in reverse order, including an empty match.
source(file.path(root_dir, "bin", "orthology_stages.R"))
permuted_groups <- tie_orthogroups
permuted_groups[, 4] <- "P4, P2, P1, P3"
permuted_groups[, 5] <- "O3, O1, O2"
stable_match <- function(groups) environment$find_ortologhi(
  canonical_orthogroups(groups), tie_preferred, tie_other,
  "ASM317683v1", "ASM207885v1", tie_final, n_amin = 3, threashold = 30)
stopifnot(identical(stable_match(permuted_groups), stable_match(tie_orthogroups)),
          identical(canonical_orthogroups(tie_orthogroups[FALSE, ]), tie_orthogroups[FALSE, ]))
second_groups <- tie_orthogroups
second_groups[, 5] <- "X1, X2, X3"
second_other <- tie_other
second_other[, 2] <- c("X1", "X2", "X3")
empty_groups <- tie_orthogroups[FALSE, ]
groups <- list(tie_orthogroups, second_groups, empty_groups)
others <- list(tie_other, second_other, tie_other)
plan <- data.frame(ordinal = 1:3, key = sprintf("genome_%06d", 1:3))
sequential <- tie_final
pairs <- list()
for (i in 1:3) {
  match_table <- function(base) environment$find_ortologhi(
    groups[[i]], tie_preferred, others[[i]], "ASM317683v1", "ASM207885v1",
    base, n_amin = 3, threashold = 30)
  sequential <- match_table(sequential)
  independent <- match_table(tie_final)
  pairs[[i]] <- list(key = plan$key[i], ordinal = i, loci = tie_final[, 2],
                    mapping = independent[, seq_len(ncol(independent)) > 4L, drop = FALSE])
}
stopifnot(identical(merge_orthology_mappings(tie_final, plan, rev(pairs)), sequential))
message("PASS: empty/tied groups and independent versus sequential orthology mapping")
