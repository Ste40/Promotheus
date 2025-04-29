### Installing and loading libraries ###

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

library(devtools)

if (!requireNamespace("rentrez", quietly = TRUE)) {
  devtools::install_github("ropensci/rentrez")
}

library(rentrez)

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

library(openxlsx)

# Load optparse for command line arguments
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}

library(optparse)

### Define command line options ###

option_list <- list(
  make_option(c("-t", "--taxID"), type = "character", default = "135461", help = "TaxID of the organism", metavar = "taxID")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Now use opt$taxID instead of hardcoding the taxID
taxID <- opt$taxID

### Function Definitions ###

get_scientific_name <- function(taxID){
  message("Retrieving the scientific name for the input taxID..")
  Taxonomy_search <- entrez_search(db = "taxonomy", term = paste0(taxID,"[UID]"))
  search_summary <- entrez_summary(db = "taxonomy", id = Taxonomy_search$ids)
  scientific_name <- search_summary$scientificname
  message(paste("The scientific name corresponding to taxID", taxID, "is", scientific_name))
  return(scientific_name)
}

get_scientific_names <- function(scientific_name){
  message(paste("Retrieving scientific names of subsp. of", scientific_name, ".."))
  scientific_name_query <- paste0(scientific_name, "[SBTR]")
  Taxonomy_search_sub <- entrez_search(db = "taxonomy", term = scientific_name_query, use_history = TRUE )
  search_summary <- entrez_summary(db = "taxonomy", web_history = Taxonomy_search_sub$web_history, retmode="xml")
  variable_names= c("Status", "Rank", "Division", "ScientificName", "CommonName", "TaxId", "AkaTaxId", "Genus", "Species", "Subsp", "ModificationDate", "GenBankDivision")
  summary_df <- as.data.frame(t((extract_from_esummary(search_summary, variable_names, simplify = TRUE))))
  scientific_names <- unname(unlist(summary_df$ScientificName))
  message(paste("Based on your input taxID, I found ", length(scientific_names), " scientific name(s) on NCBI taxonomy:\n", paste(scientific_names, collapse = "\n"), sep = ""))
  return(scientific_names)
}

# Rest of your functions go here...

# Final execution block
scientific_name <- get_scientific_name(taxID)
scientific_names <- get_scientific_names(scientific_name)
reference_genome <- get_assembly_version(scientific_name)
GEO_searches <- get_GEO_experiments_raw(scientific_names)
final_raw_db <- get_GEO_experiments_filtered(GEO_searches, scientific_names)
get_GEO_experiments_finalize(final_raw_db, scientific_names, scientific_name, reference_genome)
