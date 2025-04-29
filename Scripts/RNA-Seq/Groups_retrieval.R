
# Retrieving experimental conditions for each GSE
  
Groups_retrieval = function(GSE){ 
    gse_to_gsm_cmd= paste('pysradb gse-to-gsm',GSE, '--saveto gsm.txt')
    system(gse_to_gsm_cmd)
    gsmtable= fread("gsm.txt")
    gsm= rev(as.vector(gsmtable$experiment_alias))

    cond= list()
    
    
# Retrieving experimental conditions for each GSM
    
    for (i in 1:length(gsm)){
      
        gsmsingolo <- tryCatch(suppressMessages(getGEO(GEO= gsm[i], destdir=getwd())))

        gsm_to_srr_cmd= paste('pysradb gsm-to-srr',gsm[i], '--saveto gsm_srr.txt')
        system(gsm_to_srr_cmd)
        gsmtable= fread("gsm_srr.txt") ####continua qui, crea tabella con gsm, srr, cond e numero di appartenenza dei gruppi
        
        if(inherits(gsmsingolo, "try-error")){
            
          gsmsingolo= "error"
          
          }

      cond[[i]]= paste(c(gsmtable$experiment_alias, paste(c(gsmsingolo@header$title, gsmsingolo@header$characteristics_ch1), collapse=''), gsmtable$run_accession))
      file.remove(list.files(path=getwd(), "*.soft$"))
      
      cond[[i]]= str_replace_all(cond[[i]], "\\+", " plus ")
      cond[[i]]= str_replace_all(cond[[i]], "\\-", " minus ")
      
    }
    
    cond_dist = sapply(cond,function(x) x[2])

# Computing distances between strings to cluster experimental conditions within the same group
      
    if(length(cond_dist)>1){

      distances= adist(cond_dist)
      rownames(distances) <- colnames(distances) <- cond_dist
      threshold<- min(distances[distances > 0])+1 
      
      hc <- hclust(as.dist(distances))
      groups<- data.frame(cond_dist,cutree(hc,h=threshold)) 
      colnames(groups)<- c("Condition", "Group")

    } else { groups= data.frame(cond_dist[1], '1')
            colnames(groups)<- c("Condition", "Group")
      }

      condition_tab <- data.frame(
      GSM = character(),
      Group = character(),
      SRR = character(),
      GroupNumber = character(),
      stringsAsFactors = FALSE
      )

      for (i in 1:nrow(groups)) {
        row <- groups[i, ]
        match_condition <- row$Condition
        match_indices <- which(sapply(cond, function(x) x[2] == match_condition))
        
        for (index in match_indices) {
          to_bind <- as.data.frame(t(c(cond[[index]], row$Group)))
          colnames(to_bind) <- colnames(condition_tab)
          condition_tab[index, ] <- to_bind
        }
      }
    groups <- condition_tab
    return(groups)
  
}
  