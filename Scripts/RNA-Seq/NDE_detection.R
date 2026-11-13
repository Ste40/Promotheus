# Align labels to count columns using run IDs; never infer identity from position.
match_count_metadata <- function(files, accessions, conditions, groups) {
  ids <- sub("^([SED]RR[0-9]+)(?:_|\\.).*$", "\\1", basename(files), perl = TRUE)
  if (is.null(accessions) || length(accessions) != length(conditions) ||
      length(accessions) != length(groups) || anyNA(accessions) || anyDuplicated(accessions) ||
      anyDuplicated(ids) || !setequal(ids, accessions))
    stop("Count files and RNA-Seq run metadata do not match uniquely")
  order <- match(ids, accessions)
  list(ids = ids, conditions = conditions[order], groups = groups[order])
}

expression_statistics <- function(expression_matrix) {
  expression_matrix <- as.matrix(expression_matrix)
  storage.mode(expression_matrix) <- "double"
  means <- rowMeans(expression_matrix)
  standard_deviations <- apply(expression_matrix, 1, stats::sd)
  coefficients_of_variation <- standard_deviations / means
  coefficients_of_variation[is.na(coefficients_of_variation)] <- 0
  list(means = means, sd = standard_deviations, cv = coefficients_of_variation)
}

kruskal_pvalues <- function(expression_matrix, groups) {
  expression_matrix <- as.matrix(expression_matrix)
  storage.mode(expression_matrix) <- "double"
  if (ncol(expression_matrix) != length(groups)) {
    stop("The number of RNA-Seq samples does not match the number of group assignments")
  }
  if (length(unique(groups)) < 2L) {
    stop("At least two RNA-Seq groups are required for the Kruskal-Wallis test")
  }
  apply(expression_matrix, 1, function(values) {
    stats::kruskal.test(values, groups)$p.value
  })
}

expression_percentiles <- function(expression_matrix) {
  expression_matrix <- as.matrix(expression_matrix)
  storage.mode(expression_matrix) <- "double"
  percentiles <- apply(expression_matrix, 2, function(values) {
    rank(values, ties.method = "max") / length(values)
  })
  if (is.null(dim(percentiles))) {
    percentiles <- matrix(percentiles, ncol = 1)
  }
  percentiles
}

add_operon_status_column <- function(finaltab, status) {
  if (!is.data.frame(finaltab) || nrow(finaltab) == 0L) {
    stop("Cannot annotate operon status on an empty RNA-Seq result table")
  }

  required_columns <- c("GeneID", "locus_tag", "sd")
  missing_columns <- setdiff(required_columns, colnames(finaltab))
  if (length(missing_columns) > 0L) {
    stop("Cannot annotate operon status: RNA-Seq result table is missing columns: ", paste(missing_columns,
      collapse = ", "))
  }

  finaltab$n.genes <- rep(status, nrow(finaltab))
  finaltab
}

is_valid_rnaseq_full_tab <- function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }

  full_tab <- tryCatch(read.table(path, header = TRUE, sep = "\t", nrows = 1L), error = function(e) NULL)
  if (is.null(full_tab) || nrow(full_tab) == 0L) {
    return(FALSE)
  }

  required_columns <- c("GeneID", "locus_tag", "means", "sd", "diff_signif")
  all(required_columns %in% colnames(full_tab))
}

NDE_detection <- function(GSE, GSEdirectory, strandness, threads, m.o.subpath, GSEconditions, m.o.correct,
  od_full_tabs, CORRECTION_DE, OPERONS_DETECTION, PVALUE, run_accessions = NULL) {

  suppressWarnings(dir.create(od_full_tabs))
  setwd(GSEdirectory)

  ## Reads counting
  if (length(list.files(pattern = "*counts.txt$")) != length(list.files(pattern = "*.sorted.bam$"))) {

    bamfiles <- list.files(pattern = "*.sorted.bam$", full.names = FALSE, recursive = FALSE)
    GFF <- list.files(path = m.o.subpath, pattern = ".gff.gz$", full.names = TRUE)
    readsfiles <- list.files(pattern = ".fastq.gz$")

    if (strandness[[1]][1] == "This is SingleEnd Data") {

      if (strandness[[1]][2] > 1.2) {
        # strand+ single end
        for (i in seq_along(bamfiles)) {
          countslist <- featureCounts(files = bamfiles[i], isGTFAnnotationFile = TRUE, nthreads = threads,
          annot.ext = GFF, strandSpecific = 1, GTF.featureType = "gene", GTF.attrType = "Name",
          minMQS = 20)
          write.table(x = data.frame(countslist$annotation, countslist$counts, stringsAsFactors = FALSE),
          file = paste(readsfiles[i], "counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
        }
      } else if (strandness[[1]][2] < 0.8) {
        # reversely single end
        for (i in seq_along(bamfiles)) {
          countslist <- featureCounts(files = bamfiles[i], isGTFAnnotationFile = TRUE, nthreads = threads,
          annot.ext = GFF, strandSpecific = 2, GTF.featureType = "gene", GTF.attrType = "Name",
          minMQS = 20)
          write.table(x = data.frame(countslist$annotation, countslist$counts, stringsAsFactors = FALSE),
          file = paste(readsfiles[i], "counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
        }
      } else {
        # unstranded
        for (i in seq_along(bamfiles)) {
          countslist <- featureCounts(files = bamfiles[i], isGTFAnnotationFile = TRUE, nthreads = threads,
          annot.ext = GFF, strandSpecific = 0, GTF.featureType = "gene", GTF.attrType = "Name",
          minMQS = 20)
          write.table(x = data.frame(countslist$annotation, countslist$counts, stringsAsFactors = FALSE),
          file = paste(readsfiles[i], "counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
        }
      }

    } else {
      # Paired End
      readsfiles <- gsub(".sorted.bam", "", bamfiles)
      readsfiles <- paste0(readsfiles, ".fastq.gz")

      if (strandness[[1]][2] > 1.2) {
        # strand+ paired end
        for (i in seq_along(bamfiles)) {
          countslist <- featureCounts(files = bamfiles[i], isPairedEnd = TRUE, nthreads = threads,
          isGTFAnnotationFile = TRUE, annot.ext = GFF, strandSpecific = 1, GTF.featureType = "gene",
          GTF.attrType = "Name", minMQS = 20)
          write.table(x = data.frame(countslist$annotation, countslist$counts, stringsAsFactors = FALSE),
          file = paste(readsfiles[i], "counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
        }
      } else if (strandness[[1]][2] < 0.8) {
        # reversely paired end
        for (i in seq_along(bamfiles)) {
          countslist <- featureCounts(files = bamfiles[i], isPairedEnd = TRUE, nthreads = threads,
          isGTFAnnotationFile = TRUE, annot.ext = GFF, strandSpecific = 2, GTF.featureType = "gene",
          GTF.attrType = "Name", minMQS = 20)
          write.table(x = data.frame(countslist$annotation, countslist$counts, stringsAsFactors = FALSE),
          file = paste(readsfiles[i], "counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
        }
      } else {
        # unstranded
        for (i in seq_along(bamfiles)) {
          countslist <- featureCounts(files = bamfiles[i], isPairedEnd = TRUE, nthreads = threads,
          isGTFAnnotationFile = TRUE, annot.ext = GFF, strandSpecific = 0, GTF.featureType = "gene",
          GTF.attrType = "Name", minMQS = 20)
          write.table(x = data.frame(countslist$annotation, countslist$counts, stringsAsFactors = FALSE),
          file = paste(readsfiles[i], "counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
        }
      }
    }
  }

  ## Computing statistics and NDE identification
  output_path <- file.path(od_full_tabs, paste(GSE, m.o.correct, "full_tab.txt"))
  if (!is_valid_rnaseq_full_tab(output_path)) {

    # Normalization with DESeq2
    tablefiles <- list.files(pattern = "*fastq.gz counts.txt$", full.names = FALSE, recursive = FALSE)
    if (length(tablefiles) == 0L) {
      stop("No featureCounts tables were found for ", GSE)
    }
    final_tab <- read.table(tablefiles[1], header = TRUE, sep = "\t", quote = "")

    if (length(tablefiles) > 1L) {
      for (i in seq.int(2L, length(tablefiles))) {
        count_table <- read.table(tablefiles[i], header = TRUE, sep = "\t", quote = "")
        if (anyDuplicated(count_table$GeneID) || anyDuplicated(final_tab$GeneID) ||
            !setequal(final_tab$GeneID, count_table$GeneID)) stop("Count tables have inconsistent gene IDs")
        count_table <- count_table[match(final_tab$GeneID, count_table$GeneID), , drop = FALSE]
        final_tab <- data.frame(final_tab, count_table[, ncol(count_table)])
      }
    }

    metadata <- match_count_metadata(tablefiles, run_accessions,
      str_split(GSEconditions$Experimental_condition, ",")[[1]],
      as.numeric(strsplit(GSEconditions$Groups, ",")[[1]]))
    cond <- metadata$conditions
    tabDES <- final_tab[, -(1:6)]
    colnames(tabDES) <- cond
    rownames(tabDES) <- final_tab$GeneID

    meta <- t(as.data.frame(cond))
    colnames(meta) <- cond
    meta <- t(meta)

    dds <- suppressMessages(suppressWarnings(DESeqDataSetFromMatrix(countData = tabDES, colData = meta,
      design = ~cond)))
    dds <- suppressMessages(suppressWarnings(estimateSizeFactors(dds)))
    normalized_counts <- log(counts(dds, normalized = TRUE) + 1, 2)
    normalized_counts <- data.frame(final_tab[, 1:6], normalized_counts)

    my_data <- normalized_counts
    expression_matrix <- data.matrix(my_data[, 7:ncol(my_data), drop = FALSE])
    groups <- metadata$groups

    ## Computing gene expressions: CV, SD, means
    statistics <- expression_statistics(expression_matrix)
    means <- statistics$means
    sd <- statistics$sd
    cv <- statistics$cv
    condtabexpr <- data.frame(normalized_counts[, -c(2:6)], means, sd, cv)

    ## Kruskal-Wallis test, evaluated row-wise without rebuilding data frames.
    kruskalp <- kruskal_pvalues(expression_matrix, groups)

    ## Identifying NDE genes within a GSE
    if (CORRECTION_DE %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) {
      padjusted <- p.adjust(kruskalp, method = CORRECTION_DE, n = length(kruskalp))
    } else {
      warning("The input correction method is not valid. Using 'none' by default.")
      padjusted <- p.adjust(kruskalp, method = "none", n = length(kruskalp))
    }

    diff_signif <- ifelse(is.na(padjusted), "ZERO_COUNTS", ifelse(padjusted < PVALUE, "YES", "NO"))

    kruskalp <- data.frame(cbind(padjusted, diff_signif))
    rownames(kruskalp) <- my_data$GeneID
    difftab <- data.frame(condtabexpr, diff_signif)

    ## Computing percentiles for each expression value in each experiment
    mytab <- difftab[, seq_len(length(groups) + 1L)]
    percentile_matrix <- expression_percentiles(expression_matrix)
    percentdf <- data.frame(mytab$GeneID, percentile_matrix, check.names = FALSE)

    ## Computing percentiles mean
    permane <- rowMeans(percentile_matrix)
    genepercentiles <- data.frame(mytab$GeneID, permane)
    colnames(genepercentiles) <- c("Genes", "Percentiles")
    colnames(percentdf) <- c("GeneID", seq_along(groups))

    ## Retrieving locus tag from GFF file
    GFFtab <- readGFF(list.files(path = m.o.subpath, pattern = "*.gff.gz$", full.names = TRUE), columns = "type",
      tags = c("Name", "locus_tag"))
    GFFtab <- GFFtab[GFFtab$type == "gene", -1]
    colnames(GFFtab) <- c("GeneID", "locus_tag")

    ## Building final table
    finaltab <- cbind(difftab, genepercentiles$Percentiles, percentdf[, -1])
    colnames(finaltab)[colnames(finaltab) == "genepercentiles$Percentiles"] <- "Percentiles_mean"
    finaltab <- merge(finaltab, GFFtab, by = "GeneID")
    finaltab <- distinct(finaltab, GeneID, .keep_all = TRUE)

    ## Finding common genes between expression table (finaltab) and operon table.
    setwd(m.o.subpath)
    microorganism <- basename(m.o.subpath)
    parent_dir <- dirname(GSEdirectory)
    mode_file <- file.path(parent_dir, paste0(microorganism, "_operons_mode.txt"))

    if (!isTRUE(OPERONS_DETECTION)) {
      finaltab <- add_operon_status_column(finaltab, "OPERON_ANALYSIS_DISABLED")
    } else if (file.exists(mode_file)) {
      cat("Using operons_mode file:", mode_file, "\n")
      operons_mode <- read.table(mode_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      operons_mode$Genes <- vapply(seq_len(nrow(operons_mode)), function(i) {
        genes <- strsplit(operons_mode$Genes[i], ",\\s*", perl=TRUE)[[1]]
        if (operons_mode$Strand[i] == "+") genes[1] else genes[length(genes)]
      }, character(1))

      if (file.exists(file.path(GSEdirectory, "_operons.txt"))) {
        operonstab <- read.table(file.path(GSEdirectory, "_operons.txt"), header = TRUE, sep = "\t",
          stringsAsFactors = FALSE)
        colnames(operonstab)[ncol(operonstab)] <- "GeneID"
          operonstab$GeneID <- as.character(operonstab$GeneID)
        operonsfull <- unlist(strsplit(operonstab$GeneID, ", "))

        for (i in seq_len(nrow(operonstab))) {
          if (operonstab$Strand[i] == "+") {
          operonstab$GeneID[i] <- strsplit(operonstab$GeneID, ", ")[[i]][1]
          } else {
          gene_list <- strsplit(operonstab$GeneID, ", ")[[i]]
          operonstab$GeneID[i] <- gene_list[length(gene_list)]
          }
        }
        "%ni%" <- Negate(`%in%`)
        notfirstgenes <- operonsfull[operonsfull %ni% operonstab$GeneID]

        n.genes <- finaltab$cv * 0
        finaltab <- cbind(finaltab, n.genes)

        for (i in seq_len(nrow(finaltab))) {
          if (finaltab$GeneID[i] %in% operonstab$GeneID) {
          if (finaltab$GeneID[i] %in% operons_mode$Genes) {
            suppressWarnings(finaltab$n.genes[i] <- paste("first gene in operon (mode), n. of genes:",
            operons_mode[operons_mode$Genes == finaltab$GeneID[i], "Number.of.Genes"]))
          } else {
            suppressWarnings(finaltab$n.genes[i] <- paste("first gene in operon, n. of genes:",
            operonstab[operonstab$GeneID == finaltab$GeneID[i], "Number.of.Genes"]))
          }
          } else if (finaltab$GeneID[i] %in% notfirstgenes) {
          finaltab$n.genes[i] <- "not first gene in operon"
          } else {
          finaltab$n.genes[i] <- "gene not in operon"
          }
        }
      } else {
        finaltab <- add_operon_status_column(finaltab, "WAITING_FOR_OPERONS_TABLE")
      }
    } else {
      if (file.exists(file.path(GSEdirectory, "_operons.txt"))) {
        operonstab <- read.table(file.path(GSEdirectory, "_operons.txt"), header = TRUE, sep = "\t",
          stringsAsFactors = FALSE)
        colnames(operonstab)[ncol(operonstab)] <- "GeneID"
          operonstab$GeneID <- as.character(operonstab$GeneID)
        operonsfull <- unlist(strsplit(operonstab$GeneID, ", "))

        for (i in seq_len(nrow(operonstab))) {
          if (operonstab$Strand[i] == "+") {
          operonstab$GeneID[i] <- strsplit(operonstab$GeneID, ", ")[[i]][1]
          } else {
          gene_list <- strsplit(operonstab$GeneID, ", ")[[i]]
          operonstab$GeneID[i] <- gene_list[length(gene_list)]
          }
        }
        "%ni%" <- Negate(`%in%`)
        notfirstgenes <- operonsfull[operonsfull %ni% operonstab$GeneID]

        finaltab <- add_operon_status_column(finaltab, "WAITING_FOR_OPERONS_TABLE")
      } else {
        finaltab <- add_operon_status_column(finaltab, "WAITING_FOR_OPERONS_TABLE")
      }
    }

    finaltab <- arrange(finaltab, sd)
    write.table(finaltab, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
    token <- TRUE
    message(paste("NDE genes identified for", GSE))

  } else {
    message(paste("NDE genes table for", GSE, "already locally stored"))
    token <- TRUE
  }

  return(token)
}
