################################################################################
########################### INTER EXPERIMENT ANALYSIS ##########################
################################################################################

# Merge microarray and RNA-Seq files.
read_expression_full_tab <- function(path) {
  full_tab <- read.table(path, header = TRUE, sep = "\t")
  required_columns <- c("locus_tag", "means", "diff_signif")
  missing_columns <- setdiff(required_columns, colnames(full_tab))

  if (length(missing_columns) > 0L) {
    stop("Malformed expression full_tab '", path, "': missing columns: ", paste(missing_columns,
      collapse = ", "))
  }
  if (nrow(full_tab) == 0L) {
    stop("Malformed expression full_tab '", path, "': the table contains no genes")
  }

  full_tab
}

expression_sample_columns <- function(full_tab) {
  # MA has Geneid/locus_tag before the samples. RNA-Seq has only GeneID
  # there and appends locus_tag later, so a fixed start at column 3 loses
  # its first sample. Statistics begin at the named means column.
  mean_column <- match("means", names(full_tab))
  if (is.na(mean_column)) stop("Expression table has no means column")
  candidates <- seq_len(mean_column - 1L)
  samples <- candidates[!names(full_tab)[candidates] %in% c("Geneid", "GeneID", "locus_tag")]
  if (length(samples) == 0L || !all(vapply(full_tab[samples], is.numeric, logical(1)))) {
    stop("Expression table has no valid numeric sample columns before means")
  }
  samples
}

adjust_batch_effect <- function(expression_matrix, batches, combat = ComBat) {
  if (length(batches) != ncol(expression_matrix)) {
    stop("Expression columns and batch labels have different lengths")
  }
  if (length(unique(batches)) <= 1L) {
    return(expression_matrix)
  }

  combat(
    dat = expression_matrix,
    batch = batches,
    mod = NULL,
    par.prior = TRUE,
    prior.plots = FALSE,
    mean.only = TRUE
  )
}

unique_matrix <- function(tmp_MA_path, path_rna, tmp_genome_path, operons_detection) {

  geni_ortologhi <- read.table("/output/orthologous_genes.txt")
  # MA
  file_MA <- dir(tmp_MA_path, pattern = "_full_tab")
  gene_com <- c()
  tab <- list()
  o <- 1
  GSEcode <- c()
  if (length(file_MA) > 0) {
    for (i in seq_along(file_MA)) {
      tab[[i]] <- read_expression_full_tab(file.path(tmp_MA_path, file_MA[i]))
      if (i == 1) {
        gene_com <- tab[[i]]$locus_tag
      } else {
        gene_com <- intersect(tab[[i]]$locus_tag, gene_com)
      }
    }
    o <- length(file_MA) + 1
    GSEcode <- gsub("_full_tab.txt", "", file_MA)
  }

  # RNAseq
  file_RNAseq <- dir(path_rna, pattern = "full_tab")
  if (length(file_RNAseq) > 0) {
    for (i in seq_along(file_RNAseq)) {
      tab[[o]] <- read_expression_full_tab(file.path(path_rna, file_RNAseq[i]))
      locus_tag_unici <- tab[[o]]$locus_tag

      ind_col <- NA_integer_
      for (x in seq_along(locus_tag_unici)) {
        ind_col <- which(locus_tag_unici[x] == geni_ortologhi, arr.ind = TRUE)
        if (length(ind_col) == 2) {
          ind_col <- as.numeric(ind_col[2])
          break
        } else if (length(ind_col) > 2) {
          ind_col <- as.numeric(ind_col[1, 2])
          break
        }
      }
      if (is.na(ind_col)) {
        stop("No RNA-Seq locus tags were found in the orthology table")
      }
      riga <- c()
      for (x in seq_along(locus_tag_unici)) {
        ind <- which(locus_tag_unici[x] == geni_ortologhi[, ind_col], arr.ind = TRUE)
        if (length(ind) > 0) {
          riga[x] <- ind[1]
        } else {
          riga[x] <- NA
        }
      }
      ind <- which(is.na(riga))
      if (length(ind) > 0) {
        riga <- riga[-ind]
        tab[[o]] <- tab[[o]][-ind, ]
        locus_tag_unici <- locus_tag_unici[-ind]
      }
      tab[[o]]$locus_tag <- geni_ortologhi[riga, 2]

      gene_com <- intersect(tab[[o]]$locus_tag, gene_com)
      o <- o + 1
    }
    GSEcode <- c(GSEcode, gsub("full_tab.txt", "", file_RNAseq))
  }


  if (isTRUE(operons_detection)) {
    directories <- dir("/output/RNA_Seq")
    ind <- which(directories == "results")
    if (length(ind) > 0) {
      directories <- directories[-ind]
    }
    ind <- which(directories == "Experimental_Groups.xlsx")
    if (length(ind) > 0) {
      directories <- directories[-ind]
    }
    for (i in seq_along(directories)) {
      GSE_fulltab <- dir("/output/RNA_Seq/results")
      ind <- which(regexpr(directories[i], GSE_fulltab) > 0)
      if (length(ind) > 1) {
        ind <- ind[1]
      }
      GSE_fulltab <- read.table(file.path("/output/RNA_Seq/results", GSE_fulltab[ind]), sep = "\t",
        header = TRUE)
      locus_fulltab <- GSE_fulltab$locus_tag
      tag <- GSE_fulltab$n.genes

      ind <- which(tag == "not first gene in operon")
      if (length(ind) > 0) {
        locus_fulltab <- locus_fulltab[-ind]
      }

      ind_col <- NA_integer_
      for (x in seq_along(locus_fulltab)) {
        ind_col <- which(locus_fulltab[x] == geni_ortologhi, arr.ind = TRUE)
        if (length(ind_col) == 2) {
          ind_col <- as.numeric(ind_col[2])
          break
        } else if (length(ind_col) > 2) {
          ind_col <- as.numeric(ind_col[1, 2])
          break
        }
      }
      if (is.na(ind_col)) {
        stop("No operon locus tags were found in the orthology table")
      }

      riga <- c()
      for (x in seq_along(locus_fulltab)) {
        ind <- which(locus_fulltab[x] == geni_ortologhi[, ind_col], arr.ind = TRUE)
        if (length(ind) > 0) {
          riga[x] <- ind[1]
        } else {
          riga[x] <- NA
        }
      }

      ind <- which(is.na(riga))
      if (length(ind) > 0) {
        riga <- riga[-ind]
      }
      locus_fulltab <- geni_ortologhi[riga, 2]
      if (length(locus_fulltab) > 0) {
        gene_com <- intersect(locus_fulltab, gene_com)
        if (length(gene_com) == 0) {
          message("ERROR: all genes common to all experiments are organized in operons, therefore no results can be generated.")
        }
      }
    }
  }

  if (length(gene_com) > 0) {
    # Union expression
    ESP <- c()
    n_GSM <- c()
    DE <- rep("NO", length(gene_com))
    for (i in seq_along(tab)) {
      ind <- match(gene_com, tab[[i]]$locus_tag)
      if (anyNA(ind)) {
        stop("Internal error while aligning shared locus tags across expression tables")
      }
      sample_columns <- expression_sample_columns(tab[[i]])
      if (length(ESP) == 0) {
        ESP <- as.matrix(tab[[i]][ind, sample_columns, drop = FALSE])
      } else {
        ESP <- cbind(ESP, as.matrix(tab[[i]][ind, sample_columns, drop = FALSE]))
      }
      n_GSM[i] <- length(sample_columns)
      # nolint: trailing_whitespace_linter.  DE
      inde_DE <- which(colnames(tab[[i]]) == "diff_signif")
      ind_YES <- which(tab[[i]][ind, inde_DE] == "YES")
      if (length(ind_YES) > 0) {
        DE[ind_YES] <- "YES"
      }
    }

    gene_name <- c()
    for (i in seq_along(gene_com)) {
      ind <- which(gene_com[i] == geni_ortologhi, arr.ind = TRUE)
      if (length(ind) == 2) {
        gene_name[i] <- geni_ortologhi[ind[1], 1]
      } else if (length(ind) > 2) {
        gene_name[i] <- geni_ortologhi[ind[1, 1], 1]
      }

    }

    batches <- c(rep(GSEcode, n_GSM))
    colnames(ESP) <- batches
    final_tab <- cbind(Geneid = gene_name, locus_comuni = gene_com, ESP, DE)
    output <- list(final_tab, batches)
  } else {
    stop("No locus tags are shared by all successfully generated expression full_tab files")
  }

  return(output)
}

# Function to create final table
create_final_table <- function(espressione_media, tabella_completa) {

  gene_name <- tabella_completa[, 1]
  locus_tag <- tabella_completa[, 2]
  DE <- tabella_completa[, ncol(tabella_completa)]

  elenco_GSE <- colnames(espressione_media)

  expression_matrix <- as.matrix(espressione_media)
  storage.mode(expression_matrix) <- "double"
  mean_psnorm <- rowMeans(expression_matrix)
  sd_psnorm <- apply(expression_matrix, 1, sd)
  # A CV requires a positive mean and at least two finite observations.
  if (ncol(expression_matrix) < 2L || any(!is.finite(expression_matrix)) ||
      any(!is.finite(mean_psnorm)) || any(mean_psnorm <= 0)) {
    stop("Cannot compute CV: expression values must be finite, with at least two samples and a positive mean for every gene")
  }
  cv_psnorm <- sd_psnorm / mean_psnorm
  if (any(!is.finite(cv_psnorm))) stop("Cannot compute CV: non-finite result")

  # percentiles
  percentile <- apply(expression_matrix, 2, function(values) {
    rank(values) / nrow(expression_matrix) * 100
  })
  if (is.null(dim(percentile))) {
    percentile <- matrix(percentile, ncol = 1L)
  }
  percentile_medio <- rowMeans(percentile)
  percentile_sd <- apply(percentile, 1, sd)

  # final table
  finaltab <- cbind(gene_name, locus_tag, espressione_media, mean_psnorm, cv_psnorm, percentile, percentile_medio,
    percentile_sd, DE)
  etichette_espressione <- paste("Expression_", elenco_GSE, sep = "")
  etichette_percentile <- paste("Percentile_", elenco_GSE, sep = "")
  colnames(finaltab) <- c("Gene_name", "GeneID", etichette_espressione, "Mean_Expression", "CV", etichette_percentile,
    "Mean_Percentile", "Percentile_SD", "DE")
  finaltab <- data.frame(finaltab)

  return(finaltab)
}

# Remove temporary files
last_operations <- function(tmp_MA_path, path_rna, path_fulltab_res, path_output) {

  keep_cache <- tolower(Sys.getenv("PROMOTHEUS_KEEP_CACHE", "true")) == "true"
  files <- list.files(tmp_MA_path, pattern = "full_tab", full.names = TRUE)
  path_dest <- path_fulltab_res
  destination_paths <- file.path(path_dest, basename(files))
  if (keep_cache) {
    file.copy(files, destination_paths, overwrite = TRUE)
  } else {
    file.rename(files, destination_paths)
  }

  files <- list.files(path_rna, pattern = "full_tab", full.names = TRUE)
  path_dest <- path_fulltab_res
  destination_paths <- file.path(path_dest, gsub(" ", "_", basename(files)))
  file.copy(files, destination_paths)

  if (!keep_cache) {
    unlink(file.path(path_output, "tmp"), recursive = TRUE, force = TRUE)
  }
}

# Build current NCBI RefSeq assembly file URLs. Handles assembly-summary paths with or without a
# trailing slash.
refseq_file_location <- function(ftp_path, suffix) {
  if (length(ftp_path) != 1 || is.na(ftp_path) || !nzchar(ftp_path) || ftp_path == "na") {
    stop("RefSeq assembly has no usable ftp_path")
  }

  ftp_path <- sub("/+$", "", ftp_path)
  ftp_string <- basename(ftp_path)
  if (!nzchar(ftp_string)) {
    stop("Unable to determine RefSeq assembly name from ftp_path: ", ftp_path)
  }

  list(url = paste0(ftp_path, "/", ftp_string, suffix), filename = paste0(ftp_string, suffix))
}

promoter_sequence <- function(fasta, gene_start, gene_end, strand, promoter_length = 300L) {
  if (length(fasta) == 0L) {
    stop("Cannot extract a promoter from an empty FASTA sequence")
  }
  positions <- if (strand == "-") {
    seq.int(gene_end + 1L, gene_end + promoter_length)
  } else {
    seq.int(gene_start - promoter_length, gene_start - 1L)
  }
  circular_indices <- ((positions - 1L) %% length(fasta)) + 1L
  paste(fasta[circular_indices], collapse = "")
}

reference_gene_rows <- function(gff, gene_ids) {
  attribute_column <- intersect(c("attributes", "V9"), names(gff))
  if (length(attribute_column) == 0L) stop("GFF has no attributes column")
  attributes <- strsplit(as.character(gff[[attribute_column[1L]]]), ";", fixed = TRUE)
  loci <- vapply(attributes, function(fields) {
    tag <- fields[startsWith(fields, "locus_tag=")]
    if (length(tag) == 0L) NA_character_ else sub("^locus_tag=", "", tag[1L])
  }, character(1))
  indices <- match(as.character(gene_ids), loci)
  if (anyNA(indices)) {
    stop("GeneID not found in reference annotation: ", paste(gene_ids[is.na(indices)], collapse = ", "))
  }
  indices
}

# Function to add promter sequece (300bp) to final table
add_promoter_sequence <- function(path_output, finaltab, asm_pref) {
  summary <- fread("/output/assembly_summary.txt", quote = FALSE)

  ind <- which(summary$asm_name == asm_pref)
  ftp_path <- summary$ftp_path[ind]
  gff_location <- refseq_file_location(ftp_path, "_genomic.gff.gz")
  fna_location <- refseq_file_location(ftp_path, "_genomic.fna.gz")
  gff_path <- file.path(path_output, paste0(asm_pref, ".gff"))
  fna_path <- file.path(path_output, paste0(asm_pref, ".fna"))
  if (!file.exists(gff_path)) {
    annot <- character()
    for (i in seq_len(3L)) {
      tryCatch({
        annot <- curl_download(gff_location$url, destfile = file.path(path_output, gff_location$filename),
          mode = "wb")
        if (length(annot) > 0) {
          break
        }
      }, error = function(e) {
        message("An error occurred while downloading .gff file")
        Sys.sleep(60)
      })
    }
    if (length(annot) > 0) {
      gunzip(filename = annot, destname = gff_path, overwrite = TRUE)
    }

  }

  if (!file.exists(fna_path)) {
    annot <- character()
    for (i in seq_len(3L)) {
      tryCatch({
        annot <- curl_download(fna_location$url, destfile = file.path(path_output, fna_location$filename),
          mode = "wb")
        if (length(annot) > 0) {
          break
        }
      }, error = function(e) {
        message("An error occurred while downloading .fna file")
        Sys.sleep(60)
      })
    }
    if (length(annot) > 0) {
      gunzip(filename = annot, destname = fna_path, overwrite = TRUE)
    }

  }

  if (!file.exists(gff_path) || !file.exists(fna_path)) {
    stop("Reference GFF or FASTA file is unavailable after three download attempts")
  }

  GFF <- read.gff(gff_path)
  FASTA <- read.fasta(fna_path)[[1]]
  gene_rows <- reference_gene_rows(GFF, finaltab$GeneID)

  Sequence <- character(nrow(finaltab))
  for (i in seq_len(nrow(finaltab))) {
    ind <- gene_rows[i]
    start <- GFF[ind, 4]
    end <- GFF[ind, 5]
    strand <- as.character(GFF[ind, 7])
    Sequence[i] <- promoter_sequence(FASTA, start, end, strand)
  }
  cbind(finaltab, Sequence)
}

# Function to merge microarray and RNAseq data, normalize with Combat and create final table
interGSE_analysis <- function(path_input, tmp_MA_path, path_rna, results_path, path_fulltab_res, path_output,
  path_operons, tmp_genome_path, OPERONS_DETECTION, asm_pref, multispecies) {

  # Merge microarray and RNAseq data
  output <- unique_matrix(tmp_MA_path, path_rna, tmp_genome_path, OPERONS_DETECTION)
  if (length(output) > 0) {
    tabella_completa <- output[[1]]
    batches <- output[[2]]

    ind_col <- which(regexpr("GSE", colnames(tabella_completa)) > 0)
    E <- tabella_completa[, ind_col]
    ESP <- matrix(as.numeric(E), nrow = nrow(E), ncol = ncol(E))
    colnames(ESP) <- colnames(E)
    rm(E)

    # Divide experiment according to specie
    file_input <- dir("/input")
    GSE_tab_complete <- c()
    for (i in seq_along(file_input)) {
      text <- read.xlsx(file.path("/input", file_input[i]))
      GSE_tab <- cbind(text$Accession, text$taxon)
      if (length(GSE_tab_complete) == 0) {
        GSE_tab_complete <- GSE_tab
      } else {
        GSE_tab_complete <- rbind(GSE_tab_complete, GSE_tab)
      }
    }

    if (isTRUE(multispecies)) {
      unique_species <- unique(GSE_tab_complete[, 2])
      final_ESP <- NULL

      for (i in seq_along(unique_species)) {
        ind <- which(unique_species[i] == GSE_tab_complete[, 2])
        ESP1 <- NULL
        batches1 <- character()

        for (j in seq_along(ind)) {
          ind1 <- which(regexpr(GSE_tab_complete[ind[j], 1], colnames(ESP)) > 0)
          if (length(ind1) > 0) {
            if (is.null(ESP1)) {
              ESP1 <- ESP[, ind1, drop = FALSE]
            } else {
              ESP1 <- cbind(ESP1, ESP[, ind1, drop = FALSE])
            }
            batches1 <- c(batches1, rep(GSE_tab_complete[ind[j], 1], length(ind1)))
          }
        }

        if (is.null(ESP1)) {
          next
        }

        # INTRAEXPERIMENTAL NORMALIZATION: COMBAT
        combat_edata1 <- adjust_batch_effect(ESP1, batches1)

        if (is.null(final_ESP)) {
          final_ESP <- combat_edata1
        } else {
          final_ESP <- cbind(final_ESP, combat_edata1)
        }
      }
      if (is.null(final_ESP)) {
        stop("No expression datasets could be associated with an input species")
      }
    } else {
      batches1 <- batches
      combat_edata1 <- adjust_batch_effect(ESP, batches1)
      final_ESP <- combat_edata1
    }

    # Final table
    finaltab <- create_final_table(final_ESP, tabella_completa)

    # ADD promoter sequence
    finaltab <- add_promoter_sequence(path_output, finaltab, asm_pref)

    write.xlsx(finaltab, file.path(results_path, "final_table.xlsx"))

    # copy of full tabs file and remove of tmp directory
    last_operations(tmp_MA_path, path_rna, path_fulltab_res, path_output)
  } else {
    message("######## The analysis ended with an error. ########")
  }
}
