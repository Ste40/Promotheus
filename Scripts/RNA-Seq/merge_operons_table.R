merge_operons <- function(m.o.subpath, od_full_tabs, GSEdirectory, GSE, m.o.correct) {
  # Estrae il nome del microrganismo dalla directory
  microorganism <- basename(m.o.subpath)
  parent_dir <- dirname(GSEdirectory)
  mode_file <- file.path(parent_dir, paste0(microorganism, "_operons_mode.txt"))
  
  # Verifica se il file operons_mode è stato generato
  if(file.exists(mode_file)) {
    message("Using operons_mode file: ", mode_file)
    operons_mode <- read.table(mode_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Trova tutti i full_tab relativi a questo microorganismo nella cartella od_full_tabs.
    full_tab_files <- list.files(od_full_tabs, pattern = paste0(microorganism, ".*full_tab.txt$"), full.names = TRUE)
    if(length(full_tab_files) == 0) {
      message("No full_tab files found for microorganism ", microorganism)
      return(FALSE)
    }
    
    files_updated <- 0  # Conta quanti file sono stati aggiornati
    
    # Per ogni file full_tab, aggiorna la colonna n.genes in sostituzione dei valori "WAITING_FOR_OPERONS_TABLE"
    for(ft in full_tab_files) {
      finaltab <- read.table(ft, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      
      if(!("n.genes" %in% colnames(finaltab))) {
        message("File ", ft, " does not contain column 'n.genes'. Skipping...")
        next
      }
      
      # Procede solo per le righe in cui n.genes è ancora "WAITING_FOR_OPERONS_TABLE"
      waiting_rows <- which(finaltab$n.genes == "WAITING_FOR_OPERONS_TABLE")
      
      if(length(waiting_rows) > 0) {
        # Legge il file degli operoni copiato nella directory m.o.subpath e rinominato con GSE e m.o.correct
        new_operon_file <- file.path(m.o.subpath, paste0(GSE, " ", m.o.correct, "_operons.txt"))
        if(file.exists(new_operon_file)) {
          operonstab <- read.table(new_operon_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
          colnames(operonstab)[ncol(operonstab)] <- "GeneID"
          
          # Crea un vettore con tutti i geni elencati (separati da ", ")
          operonsfull <- unlist(strsplit(operonstab$GeneID, ", "))
          
          # Per ogni riga di operonstab, in base al senso, prende il primo (se "+") o l'ultimo gene (se "-")
          for(i in 1:nrow(operonstab)) {
            if(operonstab$Strand[i] == "+") {
              operonstab$GeneID[i] <- strsplit(operonstab$GeneID, ", ")[[i]][1]
            } else {
              gene_list <- strsplit(operonstab$GeneID, ", ")[[i]]
              operonstab$GeneID[i] <- gene_list[length(gene_list)]
            }
          }
          `%ni%` <- Negate(`%in%`)
          notfirstgenes <- operonsfull[operonsfull %ni% operonstab$GeneID]
          
          # Aggiorna la colonna n.genes in finaltab per le righe in cui il valore è "WAITING_FOR_OPERONS_TABLE"
          for(i in waiting_rows) {
            gene <- finaltab$GeneID[i]
            if(gene %in% operonstab$GeneID) {
              if(gene %in% operons_mode$Genes) {
                finaltab$n.genes[i] <- paste("first gene in operon (mode), n. of genes:",
                                              operons_mode[operons_mode$Genes == gene, "Number.of.Genes"])
              } else {
                finaltab$n.genes[i] <- paste("first gene in operon, n. of genes:",
                                              operonstab[operonstab$GeneID == gene, "Number.of.Genes"])
              }
            } else if(gene %in% notfirstgenes) {
              finaltab$n.genes[i] <- "not first gene in operon"
            } else {
              finaltab$n.genes[i] <- "gene not in operon"
            }
          }
        } else {
          message("New operon file not found in ", m.o.subpath, ". Cannot update n.genes values.")
        }
        # Riscrive il file aggiornato
        write.table(finaltab, file = ft, sep = "\t", row.names = FALSE, quote = FALSE)
        message("Updated file: ", ft)
        files_updated <- files_updated + 1
      } else {
        message("No rows in file ", ft, " require update (n.genes already set).")
      }
    }
    message("Merge of operons info completed for microorganism ", microorganism)
    return(TRUE)
  } else {
    message("operons_mode file not found for microorganism ", microorganism, ". No merge performed.")
    return(FALSE)
  }
}
