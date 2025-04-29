MetaData_retrieval <- function(microorganism_table, g, xlsx){
  
  source("/scripts/Groups_retrieval.R")
  
  GSEs = unlist(str_split(microorganism_table$GSE[g], ";")) 
  Name_ID= microorganism_table$m.o.name[g]
  Condition_list<-list()
  Acclist= list()
  
  xlsx_app <- data.frame(m.o.name= character(),
                     GSE= character(), 
                     GSM= character(),
                     Experimental_condition= character(), 
                     Groups= character(),
                     Accession_codes= character())
  
  for (r in 1:length(GSEs)){ 
    
    GSE<- GSEs[r]   
    print(GSE)
    Condition_list<- Groups_retrieval(GSE)
    Condition_list$Group<-gsub(",","",Condition_list$Group)
    
    xlsx_app[nrow(xlsx_app)+1, ] <- c(Name_ID, GSE,
                              paste(Condition_list$GSM, collapse=","), 
                              paste(Condition_list$Group, collapse=","), 
                              paste(Condition_list$GroupNumber, collapse=","),
                              paste(Condition_list$SRR, collapse=','))

  }
  return(xlsx_app)
   
}