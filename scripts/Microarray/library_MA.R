# Package installation 
packages = c("BiocManager", "GEOquery", "biomartr", "Biostrings", "rtracklayer",
             "devtools", "dplyr", "ggplot2", "reshape2", "stringr", "EnvStats", 
             "openxlsx", "curl", "data.table","limma","grDevices","base","janitor",
             "S4Vectors","stats","MultNonParam","reshape","affy","affyPLM","Biobase",
             "genefilter","BiocGenerics","arrayQualityMetrics",
             "ecolicdf","AnnotationDbi","tidyverse","RColorBrewer","car","nortest",
             "stringi","sva","genbankr","seqinr","NLP","ggpubr","utils","drc",
             "nlme","nplr","schoolmath","R.utils","ape","annotate","rentrez","readxl",
             "scales","googledrive","Rsamtools","signal","readxl","marray","swfscMisc","biomaRt",
             "GOSemSim","enrichplot","org.EcK12.eg.db","clusterProfiler","downloader","oligo") #"affyQCReport","simpleaffy"

Packages_installation = function(packages){   
  
  token = setdiff(packages, rownames(installed.packages()))
  
  if(!identical(token, character(0))) {
    
    print("I'm installing the required packages, please wait..")
    
    install.packages("BiocManager")
    # BiocManager::install("GEOquery")
    # BiocManager::install("biomartr")
    # BiocManager::install("Biostrings")
    # BiocManager::install("rtracklayer")
    
    for(i in 1:length(token)){
      BiocManager::install(token[i])
    }
  
    
    
    
    # install.packages("dplyr")
    # install.packages("ggplot2")
    # install.packages("reshape2") 
    # install.packages("stringr") 
    # install.packages("EnvStats")
    # install.packages("biomartr")
    # install.packages("openxlsx")
    # install.packages("curl")
    # install.packages("data.table")
    # install.packages("limma")
    # install.packages("grDevices")
    # install.packages("base")
    # install.packages("janitor")
    # install.packages("S4Vectors")
    # install.packages("stats")
    # install.packages("MultNonParam")
    # install.packages("reshape")
    # BiocManager::install("affy")
    # install.packages("affyPLM")
    # install.packages("Biobase")
    # install.packages("genefilter")
    # install.packages("affyQCReport")
    # install.packages("BiocGenerics")
    # install.packages("arrayQualityMetrics")
    # install.packages("simpleaffy")
    # install.packages("ecolicdf")
    # install.packages("AnnotationDbi")
    # install.packages("tidyverse")
    # install.packages("RColorBrewer")
    # install.packages("car")
    # install.packages("nortest")
    # install.packages("stringi")
    # install.packages("sva")
    # install.packages("genbankr")
    # install.packages("seqinr")
    # install.packages("NLP")
    # install.packages("ggpubr")
    # install.packages("utils")
    # install.packages("drc")
    # install.packages("nlme")
    # install.packages("nplr")
    # install.packages("schoolmath")
    # install.packages("R.utils")
    # install.packages("ape")
    # install.packages("annotate")
    # install.packages("rentrez")
    # install.packages("readxl")
    # install.packages("scales")
    # install.packages("googledrive")
    # install.packages("Rsamtools")
    # install.packages("signal")
    # install.packages("readxl")
    
    
    token = setdiff(packages, rownames(installed.packages()))
    
    if(!identical(token, character(0))){Installation = FALSE}
    
    else {
      
      Installation = TRUE
      
    } 
  } 
  
  Installation = TRUE
  

  library(limma)
  library(grDevices)
  library(base)
  library(Biostrings)
  library(janitor)
  library(S4Vectors)
  library(stats)
  library(MultNonParam)
  library(EnvStats)
  library(ggplot2)
  library(reshape)
  library(reshape2)
  library(GEOquery)
  library(affy)
  library(affyPLM)
  library(Biobase)
  library(genefilter)
  #library(affyQCReport)
  library(BiocGenerics)
  library(arrayQualityMetrics)
  #library(simpleaffy)
  library(ecolicdf)
  library(AnnotationDbi)
  library(tidyverse)
  library(dplyr)
  library(RColorBrewer)
  library(car)
  library(nortest)
  library(stringi)
  library(sva)
  library(genbankr)
  library(rtracklayer)
  library(seqinr)
  library(NLP)
  library(devtools)
  library(ggpubr)
  library(utils)
  library(drc)
  library(nlme)
  library(nplr)
  library(schoolmath)
  library(data.table)
  library(biomartr)
  library(R.utils)
  library(ape)
  library(annotate)
  library(rentrez)
  library(curl)
  library(readxl)
  library(scales)
  library(googledrive)
  library(BiocManager)
  library(Rsamtools)
  library(openxlsx)
  library(signal)
  library(readxl)
  library(stringr)
  library(marray)
  library(swfscMisc)
  library(biomaRt)
  library(GOSemSim)
  library(enrichplot)
  library(org.EcK12.eg.db)
  library(clusterProfiler)
  library(downloader)
  library(oligo)

  
  
  return(Installation)
}




