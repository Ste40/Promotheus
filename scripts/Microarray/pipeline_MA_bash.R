microarray_analysis<-function(path_input,path_output,tmp_genome_path,tmp_MA_path,CORRECTION_DE,PVALUE){
  
  file_MA<-dir(file.path(path_input),pattern = "Microarrays")

  ################################################################################
  ############################### RETRIEVE #######################################
  ################################################################################
 
  message("..Retriving GSE information from GEO NCBI..")
  source("/scripts/Microarray/metadata_retrieval.R")
  metadata_retrieval(tmp_MA_path,file_MA,path_input,path_output)

  
  ################################################################################
  ########################### INTRA EXPERIMENT ANALYSIS ##########################
  ################################################################################
  if(file.exists(file.path(tmp_MA_path,"total_metadata_GSE.xlsx"))){
    source("/scripts/Microarray/intra_experiment_analysis.R")
    source("/scripts/Microarray/one_channel_function_analysis_no_affymetrix.R")
    source("/scripts/Microarray/two_channel_function_analysis_no_affymetrix.R")
    source("/scripts/Microarray/affymetrix_analysis.R")

    
    message("..Intra-experiment analysis..")
    path_output_ort<-file.path(tmp_genome_path,"output")
    intra_exp_analysis(tmp_MA_path,tmp_genome_path, CORRECTION_DE, PVALUE)
  }else{
    message(" ################## It was not possible to analyze any of the selected experiments. ########################")
  }
  
  
}
