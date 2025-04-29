MetaData_retrieval <- function(microorganism_table, g, xlsx){
  
  source("/scripts/RNA-Seq/Groups_retrieval.R")
  GSEs=unlist(str_split(microorganism_table$GSE[g], ";"))
  GSEs=paste0("GSE", GSEs) 
  Name_ID=microorganism_table$taxon[g]
  Condition_list=list()
  Acclist=list()
  xlsx_app <- data.frame(m.o.name= character(),
                     GSE= character(), 
                     Experimental_condition= character(), 
                     Groups= character(),
                     Accession_codes= character())
  
  for (r in 1:length(GSEs)){ 
    
    GSE<- GSEs[r]
    Condition_list[[r]]<- Groups_retrieval(GSE)
    Condition_list[[r]]<- append(GSE, Condition_list[[r]])
    Condition_list[[r]]<- append(Name_ID, Condition_list[[r]])

    Condition_list[[r]]$Group<-gsub(",","",Condition_list[[r]]$Group)
    
    srp_to_gse_cmd<- paste('pysradb gse-to-srp',GSE, '--saveto gse.txt')
    system(srp_to_gse_cmd)
    srp_table<- fread("gse.txt")
    srp<- srp_table[1,"study_accession"]
    
    srp_metadata<- paste('pysradb metadata',srp,'--saveto srp_metadata.txt')
    system(srp_metadata)
    metadata_table<- fread("srp_metadata.txt")

    srp_to_srr_cmd<- paste('pysradb srp-to-srr',srp,'--saveto srp.txt')
    system(srp_to_srr_cmd)
    
    full_report<- fread("srp.txt")
    Acclist[[r]]<- rev(full_report$run_accession)
    xlsx_app[nrow(xlsx_app)+1, ] <- c(unlist(Condition_list[[r]][1]), 
                              unlist(Condition_list[[r]][2]), 
                              paste(as.vector(unlist(lapply(Condition_list[r], function(x) x['Group']))), collapse=","), 
                              paste(as.vector(unlist(lapply(Condition_list[r], function(x) x['GroupNumber']))), collapse=","),
                              paste(Acclist[[r]], collapse=','))

  }
  return(xlsx_app)
}