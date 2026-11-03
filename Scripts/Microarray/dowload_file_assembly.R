# Build current NCBI RefSeq assembly file URLs. The assembly summary now returns
# ftp_path values with a trailing slash, while older snapshots did not.
refseq_file_location <- function(ftp_path, suffix) {
  if (length(ftp_path) != 1 || is.na(ftp_path) || !nzchar(ftp_path) || ftp_path == "na") {
    stop("RefSeq assembly has no usable ftp_path")
  }

  ftp_path <- sub("/+$", "", ftp_path)
  ftp_string <- basename(ftp_path)
  if (!nzchar(ftp_string)) {
    stop("Unable to determine RefSeq assembly name from ftp_path: ", ftp_path)
  }

  list(
    url = paste0(ftp_path, "/", ftp_string, suffix),
    filename = paste0(ftp_string, suffix)
  )
}

# Download a compressed GenBank assembly and extract it locally.
download <- function(ftp_path, path_dir_org, asm) {
    location <- refseq_file_location(ftp_path, "_genomic.gbff.gz")
    annot <- curl_download(
      location$url,
      destfile = file.path(path_dir_org, location$filename),
      mode = "wb"
    )
    new_name <- file.path(path_dir_org, paste0(asm, ".gb"))
    gunzip(filename = annot, destname = new_name, overwrite = TRUE)

}

# Function to read assembly_summary and to download file for each version of genome
reading_testo <- function(path_output, tmp_genome_path, asm_code, asm_pref) {

  summary_path <- file.path(path_output, "assembly_summary.txt")
  if (!file.exists(summary_path)){
      message("Downloading assembly summary from RefSeq..")
      curl_download(
        "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt",
        destfile = summary_path,
        mode = "wb"
      )
  }

  summary <- fread(summary_path, quote = FALSE)

  requested_assemblies <- unique(c(asm_code, asm_pref))
  assembly_indices <- match(requested_assemblies, summary$asm_name)
  if (anyNA(assembly_indices)) {
    stop(
      "Assembly name(s) not found in the current RefSeq summary: ",
      paste(requested_assemblies[is.na(assembly_indices)], collapse = ", ")
    )
  }

  taxid_species <- unique(summary$species_taxid[assembly_indices])
  row <- summary[summary$species_taxid %in% taxid_species,]
  row <- row[!is.na(row$ftp_path) & row$ftp_path != "na" & nzchar(row$ftp_path),]
  if (nrow(row) == 0L) {
    stop("No downloadable RefSeq assemblies were found for the requested species")
  }

  ind<-which(row$assembly_level=="Complete Genome")
  if(length(ind)>0){
    row<-row[ind,]
  }

  strain_name<-c()
  for(i in 1:nrow(row)){
    strain_name[i]<-paste(row$organism_name[i],gsub("="," ",row$infraspecific_name[i]))
    strain_name[i]<-gsub(" ","_",strain_name[i])
    strain_name[i]<-gsub("\\.","",strain_name[i])
  }
  unique_strain_name<-unique(strain_name)

  for(i in seq_along(unique_strain_name)){
    ind<-which(unique_strain_name[i]==strain_name)
    if(length(ind)>1){
      ind<-ind[1]
    }
    ftp_path<-row$ftp_path[ind]
    asm<-row$asm_name[ind]
    download(ftp_path,tmp_genome_path,asm)
  }
}
