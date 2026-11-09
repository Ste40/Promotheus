# Retain the GSM -> SRR association instead of joining independent API orders.
canonical_run_conditions <- function(tab) {
  required <- c("GSM", "SRR", "Group", "GroupNumber")
  if (!all(required %in% names(tab)) || !nrow(tab)) stop("Missing RNA-Seq run metadata")
  if (anyNA(tab[required]) || any(!grepl("^[SED]RR[0-9]+$", tab$SRR)) ||
      anyDuplicated(tab$SRR) || anyDuplicated(tab$GSM))
    stop("RNA-Seq requires unique run IDs and one run per GSM; multiple runs per sample are not supported")
  groups <- suppressWarnings(as.integer(as.character(tab$GroupNumber)))
  if (anyNA(groups) || any(groups < 1L)) stop("Invalid RNA-Seq group assignments")
  tab[order(as.character(tab$SRR), method = "radix"), , drop = FALSE]
}

MetaData_retrieval <- function(microorganism_table, g, xlsx) {
  source("/scripts/RNA-Seq/Groups_retrieval.R")
  GSEs <- paste0("GSE", unlist(str_split(microorganism_table$GSE[g], ";")))
  rows <- lapply(GSEs, function(gse) {
    tab <- canonical_run_conditions(Groups_retrieval(gse))
    data.frame(m.o.name = microorganism_table$taxon[g], GSE = gse,
      Experimental_condition = paste(gsub(",", "", tab$Group), collapse = ","),
      Groups = paste(tab$GroupNumber, collapse = ","),
      Accession_codes = paste(tab$SRR, collapse = ","), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
