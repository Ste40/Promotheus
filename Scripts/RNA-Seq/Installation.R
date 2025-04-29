
# Packages installation #
# This script contains all the packages required from the pipeline #
# if you want to add a new package, remember to update the vector packages and the installation within the function #  

packages = c(# RNA-Seq #
              "BiocManager", "Rsubread", "GEOquery", "biomaRt", "Biostrings", "DESeq2", "rtracklayer",
             "remotes", "devtools", "ProGenome", "dplyr", "ggplot2", "reshape2", "stringr", "EnvStats", 
             "openxlsx", "curl", "data.table", 
             # Microarrays # 
             "janitor", "MultNonParam", "affy", "affyPLM", "genefilter",
             "arrayQualityMetrics", "ecolicdf", "tidyverse", "RColorBrewer", "car", "sva", "genbankr", "seqinr",
             "NLP", "ggpubr", "drc", "ape", "annotate", "rentrez", "signal", "readxl" )

Packages_installation = function(packages){   

  token = setdiff(packages, rownames(installed.packages()))

  if(!identical(token, character(0))) {
    
    print("I'm installing the required packages, please wait..")

  install.packages("BiocManager")
    BiocManager::install(c("Rsubread", "GEOquery", "Biostrings", "DESeq2", 
                          "rtracklayer", "biomaRt", "affy", "affyPLM",
                          "genefilter", "arrayQualityMetrics","ecolicdf", 
                          "sva", "genbankr", "annotate"), update = FALSE, ask = FALSE)
  
  install.packages(c("ggplot2", "reshape2", "stringr", "EnvStats",
                    "openxlsx", "data.table", "biomartr", "janitor",
                    "MultNonParam", "tidyverse", "dplyr", "curl",
                    "RColorBrewer", "car", "seqinr", "NLP", "ggpubr",
                    "drc", "ape", "rentrez", "signal", "readxl"), force = TRUE, quiet = TRUE)

    
  install.packages("remotes")
    remotes::install_github("YulongNiu/ParaMisc")
  
  install.packages("devtools")
    devtools::install_github('YulongNiu/KEGGAPI')
    devtools::install_github('YulongNiu/NCBIAPI')
    devtools::install_github('YulongNiu/ProGenome')
    
  
  token = setdiff(packages, rownames(installed.packages()))
  
  if(!identical(token, character(0))){Installation = FALSE}
  
   else {
    
    Installation = TRUE
    
    } 
  } 
  
  Installation = TRUE
  return(Installation)
}



  
  
  




