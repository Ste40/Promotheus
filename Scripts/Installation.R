
# Packages installation #
# This script contains all the packages required from the pipeline #
# if you want to add a new package, remember to update the vector packages and the installation within the function #  

packages = c(# RNA-Seq #
              "BiocManager", "Rsubread", "GEOquery", "biomaRt", "Biostrings", "DESeq2", "rtracklayer",
             "remotes", "devtools", "ProGenome", "dplyr", "ggplot2", "reshape2", "stringr", "EnvStats", 
             "openxlsx", "curl", "data.table", 
             # Microarrays # 
             "janitor", "MultNonParam", "affy", "affyPLM", "genefilter", "nplr", "limma", 
             "S4Vectors", "Biobase", "BiocGenerics", "AnnotationDbi", "nortest", "stringi", 
             "arrayQualityMetrics", "ecolicdf", "tidyverse", "RColorBrewer", "car", "sva", "seqinr",
             "NLP", "ggpubr", "drc", "ape", "annotate", "rentrez", "signal", "readxl", "nlme", "schoolmath",
             "R.utils", "scales", "googledrive", "Rsamtools", "marray", "swfscMisc", "GOSemSim", "enrichplot", 
             "org.EcK12.eg.db", "clusterProfiler", "downloader", "oligo", "genbankr")


Packages_installation = function(packages){   

  token = setdiff(packages, rownames(installed.packages()))
  
  if(!identical(token, character(0))) {

    message(paste("The following packages are missing:", token))
    
    print("I'm installing the required packages, please wait..")

  install.packages("BiocManager")
    BiocManager::install(c("Rsubread", "GEOquery", "Biostrings", "DESeq2", 
                          "rtracklayer", "biomaRt", "affy", "affyPLM",
                          "genefilter", "arrayQualityMetrics","ecolicdf", 
                          "sva", "annotate", "limma", "S4Vectors", 
                          "Biobase", "BiocGenerics", "AnnotationDbi", "Rsamtools", 
                          "marray", "GOSemSim", "enrichplot", "org.EcK12.eg.db", 
                          "clusterProfiler"), update = FALSE, ask = FALSE)
    
    BiocManager::install("oligo", configure.args="--disable-threading", force = TRUE)
  
    install.packages(c("ggplot2", "reshape2", "stringr", "EnvStats",
                    "openxlsx", "data.table", "biomartr", "janitor",
                    "MultNonParam", "tidyverse", "dplyr", "curl",
                    "RColorBrewer", "car", "seqinr", "NLP", "ggpubr",
                    "drc", "ape", "rentrez", "signal", "readxl", "nplr",
                    "nortest", "stringi", "nlme", "schoolmath", "R.utils",
                    "scales", "googledrive", "swfscMisc", "downloader"), force = TRUE, quiet = TRUE)

    
  install.packages("remotes")
    remotes::install_github("YulongNiu/ParaMisc")
    remotes::install_github("gmbecker/genbankr")
  
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

message("Checking if the required packages are installed ")

if(Packages_installation(packages)==FALSE){
  stop("Oooops! Something went wrong during package installation")
} else {message ("Installation successful!")}


  
  
  




