Operons_detection <- function(m.o.subpath, threads, memory, Acclist, GSEdirectory, GSEconditions, temp_threshold, gse_count) {
  
  # Estrae il nome del microrganismo dalla directory m.o.subpath
  microorganism <- basename(m.o.subpath)
  
  # Ottiene il numero di GSE per il microrganismo dal data.frame gse_count
  gse_for_mo <- as.numeric(gse_count$GSE[str_replace_all(gse_count$m.o.name, "[[:punct:]]", " ") == microorganism])
  if(length(gse_for_mo) == 0) {
    warning("Microorganism ", microorganism, " not found in gse_count. Using user-provided operon_threshold.")
  } else if(temp_threshold > gse_for_mo) {
    warning("User-defined operon_threshold (", temp_threshold, ") is greater than available GSE (", gse_for_mo, ") for microorganism ", microorganism, ". Reassigning operon_threshold to ", gse_for_mo, ".")
    temp_threshold <- gse_for_mo
  }
  
    # Se il file degli operoni non esiste ancora nella directory corrente, procedi con la rilevazione
  if(identical(list.files(path = GSEdirectory, pattern = "*_operons.txt$"), character(0))) {
    
    Rockhopperpath <- "/scripts/RNA-Seq/Rockhopper.jar"
    
    readsfiles <- list.files(pattern = "*.fastq.gz$", full.names = FALSE, recursive = FALSE, path = GSEdirectory)
    reads <- c()
    
    if(length(readsfiles) != length(Acclist)) {
      Groups <- as.data.frame(str_split(GSEconditions, ","))
      k <- 1
      for(i in seq(1, length(readsfiles), by = 2)) {
        reads[k] <- paste(readsfiles[i], readsfiles[i+1], sep = "%")
        k <- k + 1
      }
      
      dfgroups <- data.frame(reads, Groups)
      groupscounter <- 1
      groupsvector <- c()
      groupsvector[1] <- dfgroups$reads[1]
      for(i in 2:nrow(dfgroups)) {
        if(dfgroups[i, 2] == groupscounter) {
          groupsvector[groupscounter] <- paste0(groupsvector[groupscounter], ",", dfgroups[i, 1])
        } else {
          groupscounter <- groupscounter + 1
          groupsvector[groupscounter] <- dfgroups[i, 1]
        }
      }
      groupsvector <- paste(groupsvector, collapse = ' ')
      
    } else {
      Groups <- as.data.frame(str_split(GSEconditions, ","))
      dfgroups <- data.frame(readsfiles, Groups)
      groupscounter <- 1
      groupsvector <- c()
      groupsvector[1] <- readsfiles[1]
      for(i in 2:nrow(dfgroups)) {
        if(dfgroups[i, 2] == groupscounter) {
          groupsvector[groupscounter] <- paste0(groupsvector[groupscounter], ",", dfgroups[i, 1])
        } else {
          groupscounter <- groupscounter + 1
          groupsvector[groupscounter] <- dfgroups[i, 1]
        }
      }
      groupsvector <- paste(groupsvector, collapse = ' ')
    }
    
    # Esecuzione di Rockhopper
    setwd(GSEdirectory)
    genomespath <- gsub(" ", "\\\\ ", m.o.subpath)
    
    if (identical(list.files(path = m.o.subpath, pattern = "*.fna$"), character(0))) {
      string <- paste("gunzip", paste0(genomespath, "/", list.files(pattern = "*.fna.gz$", path = m.o.subpath)))
      bash <- paste('bash -c', shQuote(string))
      system(bash)
    }
    
    outpath <- GSEdirectory
    outpath <- gsub(" ", "\\\\ ", GSEdirectory)
    Rockhopperpath <- gsub(" ", "\\\\ ", Rockhopperpath)
    
    string <- paste("java", paste0("-Xmx", memory, "m"), "-cp", Rockhopperpath, "Rockhopper -o", outpath, "-g", genomespath, groupsvector)
    bash <- paste('bash -c', shQuote(string))
    system(bash)
    
  }

  parent_dir <- dirname(GSEdirectory)
  gse_dirs <- list.dirs(path = parent_dir, full.names = TRUE, recursive = FALSE)
  detected_count <- 0
  for(gse in gse_dirs) {
    if(length(list.files(path = gse, pattern = "*_operons.txt$")) > 0) {
      detected_count <- detected_count + 1
    }
  }
  
mode_filename <- file.path(parent_dir, paste0(microorganism, "_operons_mode.txt"))
if(file.exists(mode_filename)) {
  cat("Operon detection already performed for", temp_threshold, 
      "different GSE datasets for this microorganism. Skipping detection.\n")
  return(TRUE)
} else if(detected_count == temp_threshold) {
  # Legge tutti i file *_operons.txt presenti nelle directory dei GSE
  all_data <- NULL
  for(gse in gse_dirs) {
    op_files <- list.files(path = gse, pattern = "*_operons.txt$", full.names = TRUE)
    if(length(op_files) > 0) {
      df <- read.table(op_files[1], header = TRUE, sep = "\t", 
                       stringsAsFactors = FALSE, check.names = FALSE)
      if(is.null(all_data)) {
        all_data <- df
      } else {
        all_data <- rbind(all_data, df)
      }
    }
  }
  
  if(!is.null(all_data)) {
    library(dplyr)
    all_data <- all_data %>% 
      mutate(row_str = paste(`Start`, `Stop`, `Strand`, `Number of Genes`, Genes, sep = "\t"))
    
    mode_data <- all_data %>%
      group_by(Genes, row_str) %>%
      summarise(freq = n(), .groups = "drop") %>%
      group_by(Genes) %>%
      slice_max(freq, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    mode_data <- mode_data %>%
      mutate(Start = sapply(strsplit(row_str, "\t"), `[`, 1),
             Stop = sapply(strsplit(row_str, "\t"), `[`, 2),
             Strand = sapply(strsplit(row_str, "\t"), `[`, 3),
             `Number of Genes` = sapply(strsplit(row_str, "\t"), `[`, 4)) %>%
      select(Start, Stop, Strand, `Number of Genes`, Genes)
    
    write.table(mode_data, file = mode_filename, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Mode operons file created:", mode_filename, "\n")
  } else {
    cat("No operons data found to generate mode file.\n")
  }
}
  
  if(file.exists(file.path(GSEdirectory, "_operons.txt"))) {
      new_filename <- file.path(m.o.subpath, paste0(GSE, " ", m.o.correct, "_operons.txt"))
      file.copy(file.path(GSEdirectory, "_operons.txt"), new_filename, overwrite = TRUE)
      cat("Operons file copied to:", new_filename, "\n")
    token <- TRUE
  } else {
    token <- FALSE
  }
  
  return(token)
}
