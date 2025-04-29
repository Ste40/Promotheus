NDE_detection <- function(GSEdirectory, strandness, threads, m.o.subpath, GSEconditions, m.o.correct, od_full_tabs, CORRECTION_DE, OPERONS_DETECTION, PVALUE) {
  
  suppressWarnings(dir.create(od_full_tabs))
  setwd(GSEdirectory)
  
  ## Reads counting
  if(length(list.files(pattern="*counts.txt$")) != length(list.files(pattern="*.sorted.bam$"))) {
    
    bamfiles <- list.files(pattern="*.sorted.bam$", full.names=FALSE, recursive=FALSE)
    GFF <- list.files(path=m.o.subpath, pattern=".gff.gz$", full.names=TRUE)
    readsfiles <- list.files(pattern=".fastq.gz$")
    
    if(strandness[[1]][1] == "This is SingleEnd Data") {
      
      if(strandness[[1]][2] > 1.2) { # strand+ single end
        for(i in 1:length(bamfiles)) {
          countslist <- featureCounts(files = bamfiles[i], isGTFAnnotationFile = TRUE, nthreads = threads,
                                      annot.ext = GFF, strandSpecific = 1, GTF.featureType = "gene", 
                                      GTF.attrType = "Name", minMQS = 20)
          write.table(x = data.frame(countslist$annotation, countslist$counts, stringsAsFactors = FALSE),
                      file = paste(readsfiles[i], "counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
        }
      } else if(strandness[[1]][2] < 0.8) { # reversely single end
        for(i in 1:length(bamfiles)) {
          countslist <- featureCounts(files = bamfiles[i], isGTFAnnotationFile = TRUE, nthreads = threads,
                                      annot.ext = GFF, strandSpecific = 2, GTF.featureType = "gene", 
                                      GTF.attrType = "Name", minMQS = 20)
          write.table(x = data.frame(countslist$annotation, countslist$counts, stringsAsFactors = FALSE),
                      file = paste(readsfiles[i], "counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
        }
      } else { # unstranded
        for(i in 1:length(bamfiles)) {
          countslist <- featureCounts(files = bamfiles[i], isGTFAnnotationFile = TRUE, nthreads = threads,
                                      annot.ext = GFF, strandSpecific = 0, GTF.featureType = "gene", 
                                      GTF.attrType = "Name", minMQS = 20)
          write.table(x = data.frame(countslist$annotation, countslist$counts, stringsAsFactors = FALSE),
                      file = paste(readsfiles[i], "counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
        }
      }
      
    } else { # Paired End
      readsfiles <- gsub(".sorted.bam", "", bamfiles)
      readsfiles <- paste0(readsfiles, ".fastq.gz")
      
      if(strandness[[1]][2] > 1.2) { # strand+ paired end
        for(i in 1:length(bamfiles)) {
          countslist <- featureCounts(files = bamfiles[i], isPairedEnd = TRUE, nthreads = threads, 
                                      isGTFAnnotationFile = TRUE, annot.ext = GFF, strandSpecific = 1, 
                                      GTF.featureType = "gene", GTF.attrType = "Name", minMQS = 20)
          write.table(x = data.frame(countslist$annotation, countslist$counts, stringsAsFactors = FALSE),
                      file = paste(readsfiles[i], "counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
        }
      } else if(strandness[[1]][2] < 0.8) { # reversely paired end
        for(i in 1:length(bamfiles)) {
          countslist <- featureCounts(files = bamfiles[i], isPairedEnd = TRUE, nthreads = threads, 
                                      isGTFAnnotationFile = TRUE, annot.ext = GFF, strandSpecific = 2, 
                                      GTF.featureType = "gene", GTF.attrType = "Name", minMQS = 20)
          write.table(x = data.frame(countslist$annotation, countslist$counts, stringsAsFactors = FALSE),
                      file = paste(readsfiles[i], "counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
        }
      } else { # unstranded
        for(i in 1:length(bamfiles)) {
          countslist <- featureCounts(files = bamfiles[i], isPairedEnd = TRUE, nthreads = threads, 
                                      isGTFAnnotationFile = TRUE, annot.ext = GFF, strandSpecific = 0, 
                                      GTF.featureType = "gene", GTF.attrType = "Name", minMQS = 20)
          write.table(x = data.frame(countslist$annotation, countslist$counts, stringsAsFactors = FALSE),
                      file = paste(readsfiles[i], "counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
        }
      }
    }
  }
  
  ## Computing statistics and NDE identification
  if(identical(list.files(pattern = paste(GSE, m.o.correct, "full_tab.txt"), path = m.o.subpath), character(0))) {
    
    # Normalization with DESeq2
    tablefiles <- list.files(pattern = "*fastq.gz counts.txt$", full.names = FALSE, recursive = FALSE)
    final_tab <- read.table(tablefiles[1], header = TRUE, sep = "\t", quote = "")
    
    for(i in 2:length(tablefiles)) {
      count_table <- read.table(tablefiles[i], header = TRUE, sep = "\t", quote = "")
      final_tab <- data.frame(final_tab, count_table[, ncol(count_table)])
    }
    
    cond <- str_split(GSEconditions$Experimental_condition, ",")[[1]]
    tabDES <- final_tab[, -(1:6)]
    colnames(tabDES) <- cond
    rownames(tabDES) <- final_tab$GeneID
    
    meta <- t(as.data.frame(cond))
    colnames(meta) <- cond
    meta <- t(meta)
    
    dds <- suppressMessages(suppressWarnings(DESeqDataSetFromMatrix(countData = tabDES, colData = meta, design = ~ cond)))
    dds <- suppressMessages(suppressWarnings(estimateSizeFactors(dds)))
    normalized_counts <- log(counts(dds, normalized = TRUE) + 1, 2)
    normalized_counts <- data.frame(final_tab[, 1:6], normalized_counts)
    
    ## Computing gene expressions: CV, SD, means
    kruskalp <- c()
    sd <- c()
    means <- c()
    cv <- c()
    my_data <- normalized_counts
    
    for(i in 1:length(my_data$GeneID)) {
      singlegene <- my_data[i, c(1, 7:ncol(my_data))]
      singlegenetab <- as.data.frame(t(singlegene))
      singlegenetab <- cbind(rownames(singlegenetab), singlegenetab[1])
      singlegenetab <- singlegenetab[-1, ]
      rownames(singlegenetab) <- 1:nrow(singlegenetab)
      colnames(singlegenetab) <- c("groups", "expr")
      singlegenetab$groups <- as.numeric(strsplit(GSEconditions$Groups, ",")[[1]])
      
      means[i] <- mean(as.numeric(singlegenetab$expr))
      sd[i] <- sd(singlegenetab$expr)
      cv[i] <- sd[i] / means[i]
    }
    
    cv[is.na(cv)] <- 0
    condtabexpr <- data.frame(normalized_counts[, -c(2:6)], means, sd, cv)
    
    ## Fixing data for Kruskal Wallis test
    kruskalp <- c()
    for(i in 1:length(my_data$GeneID)) {
      singlegene <- my_data[i, c(1, 7:ncol(my_data))]
      singlegenetab <- as.data.frame(t(singlegene))
      singlegenetab <- cbind(rownames(singlegenetab), singlegenetab[1])
      singlegenetab <- singlegenetab[-1, ]
      rownames(singlegenetab) <- 1:nrow(singlegenetab)
      colnames(singlegenetab) <- c("groups", "expr")
      singlegenetab$groups <- as.numeric(strsplit(GSEconditions$Groups, ",")[[1]])
      
      kruskalp[i] <- kruskal.test(expr ~ groups, data = singlegenetab)$p.value
    }
    
    ## Identifying NDE genes within a GSE
    if(CORRECTION_DE %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) {
      padjusted <- p.adjust(kruskalp, method = CORRECTION_DE, n = length(kruskalp))
    } else {
      warning("The input correction method is not valid. Using 'none' by default.")
      padjusted <- p.adjust(kruskalp, method = "none", n = length(kruskalp))
    }
    
    diff_signif <- c()
    for(i in 1:length(kruskalp)) {
      if(is.numeric(padjusted[i]) && !is.nan(padjusted[i]) && !is.na(padjusted[i]) && padjusted[i] < PVALUE) {
        diff_signif[i] <- "YES"
      } else if(is.nan(padjusted[i]) || is.na(padjusted[i])) {
        diff_signif[i] <- "ZERO_COUNTS"
      } else {
        diff_signif[i] <- "NO"
      }
    }
    
    kruskalp <- data.frame(cbind(padjusted, diff_signif))
    rownames(kruskalp) <- my_data$GeneID
    difftab <- data.frame(condtabexpr, diff_signif)
    
    ## Computing percentiles for each expression value in each experiment
    percentiles <- c()
    mytab <- difftab[, 0:(length(strsplit(GSEconditions$Groups, ",")[[1]]) + 1)]
    percentdf <- mytab$GeneID
    for(i in 2:ncol(mytab)) {
      ordereddistrib <- sort(mytab[, i])
      for(j in 1:nrow(mytab)) {
        percentiles[j] <- ecdf(ordereddistrib)(mytab[j, i])
      }
      percentdf <- data.frame(percentdf, percentiles)
    }
    
    ## Computing percentiles mean
    permane <- c()
    for(i in 1:nrow(percentdf)) {
      permane[i] <- mean(as.numeric(percentdf[i, -1]))
    }
    genepercentiles <- data.frame(mytab$GeneID, permane)
    colnames(genepercentiles) <- c("Genes", "Percentiles")
    colnames(percentdf) <- c("GeneID", 1:length(strsplit(GSEconditions$Groups, ",")[[1]]))
    
    ## Retrieving locus tag from GFF file
    GFFtab <- readGFF(list.files(path = m.o.subpath, pattern = "*.gff.gz$", full.names = TRUE), columns = "type", tags = c("Name", "locus_tag"))
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
    
    if(file.exists(mode_file)) {
      cat("Using operons_mode file:", mode_file, "\n")
      operons_mode <- read.table(mode_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      
      if(file.exists(file.path(GSEdirectory, "_operons.txt"))) {
        operonstab <- read.table(file.path(GSEdirectory, "_operons.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        colnames(operonstab)[ncol(operonstab)] <- "GeneID"
        operonsfull <- unlist(strsplit(operonstab$GeneID, ", "))
        
        for(i in 1:nrow(operonstab)) {
          if(operonstab$Strand[i] == "+") {
            operonstab$GeneID[i] <- strsplit(operonstab$GeneID, ", ")[[i]][1]
          } else {
            gene_list <- strsplit(operonstab$GeneID, ", ")[[i]]
            operonstab$GeneID[i] <- gene_list[length(gene_list)]
          }
        }
        '%ni%' <- Negate(`%in%`)
        notfirstgenes <- operonsfull[operonsfull %ni% operonstab$GeneID]
        
        n.genes <- finaltab$cv * 0
        finaltab <- cbind(finaltab, n.genes)
        
        for(i in 1:nrow(finaltab)) {
          if(finaltab$GeneID[i] %in% operonstab$GeneID) {
            if(finaltab$GeneID[i] %in% operons_mode$Genes) {
              suppressWarnings(
                finaltab$n.genes[i] <- paste("first gene in operon (mode), n. of genes:", 
                                             operons_mode[operons_mode$Genes == finaltab$GeneID[i], "Number.of.Genes"])
              )
            } else {
              suppressWarnings(
                finaltab$n.genes[i] <- paste("first gene in operon, n. of genes:", 
                                             operonstab[operonstab$GeneID == finaltab$GeneID[i], "Number.of.Genes"])
              )
            }
          } else if(finaltab$GeneID[i] %in% notfirstgenes) {
            finaltab$n.genes[i] <- "not first gene in operon"
          } else {
            finaltab$n.genes[i] <- "gene not in operon"
          }
        }
      } else {
        n.genes <- rep("WAITING_FOR_OPERONS_TABLE", nrow(finaltab))
        finaltab <- cbind(finaltab, n.genes)
      }
    } else {
      if(file.exists(file.path(GSEdirectory, "_operons.txt"))) {
        operonstab <- read.table(file.path(GSEdirectory, "_operons.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        colnames(operonstab)[ncol(operonstab)] <- "GeneID"
        operonsfull <- unlist(strsplit(operonstab$GeneID, ", "))
        
        for(i in 1:nrow(operonstab)) {
          if(operonstab$Strand[i] == "+") {
            operonstab$GeneID[i] <- strsplit(operonstab$GeneID, ", ")[[i]][1]
          } else {
            gene_list <- strsplit(operonstab$GeneID, ", ")[[i]]
            operonstab$GeneID[i] <- gene_list[length(gene_list)]
          }
        }
        '%ni%' <- Negate(`%in%`)
        notfirstgenes <- operonsfull[operonsfull %ni% operonstab$GeneID]
        
        n.genes <- rep("WAITING_FOR_OPERONS_TABLE", nrow(finaltab))
        finaltab <- cbind(finaltab, n.genes)
      } else {
        finaltab <- data.frame(n.genes = rep("WAITING_FOR_OPERONS_TABLE", nrow(finaltab)), stringsAsFactors = FALSE)
      }
    }
    
    finaltab <- arrange(finaltab, sd)
    write.table(finaltab, file = file.path(od_full_tabs, paste(GSE, m.o.correct, "full_tab.txt")), sep = "\t", row.names = FALSE, quote = FALSE)
    token <- TRUE
    message(paste("NDE genes identified for", GSE))
    
  } else {
    message(paste("NDE genes table for", GSE, "already locally stored"))
    token <- TRUE
  }
  
  return(token)
}
