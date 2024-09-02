#' tflamingo.data_prepare
#'
#' Prepare data for tensor completion
#' @param code_path code path of tFlamingorLite
#' @param input_folder Path to the scHiC data folder.
#' @param low_res low resolutuion backbone
#' @param high_res high resolution
#' @param outputs_folder Path to the folder to store intermediate outputs.
#' @param chr_name Name of the chromosome
#' @param assembly Genome assembly version
#' @export


tflamingo.data_prepare = function(code_path, input_folder, low_res, high_res, outputs_folder,chr_name, assembly){
  library(tFlamingorLite)
  library(GenomeInfoDb)
  # first generate contact maps at low-resolution (e.g.300kb)
  tflamingo.generate_matrix(input_path=input_folder, resolution=low_res, chr_name = chr_name, opt_path=paste0(outputs_folder,'/lowres_contact_maps'),assembly=assembly)
  print("Finished generating matrix for low resolution")
  # then generate contact maps at high-resolution (e.g.30kb)
  tflamingo.generate_matrix(input_path=input_folder, resolution=high_res, chr_name = chr_name, opt_path=paste0(outputs_folder,'/highres_contact_maps'),assembly=assembly)
  print("Finished generating matrix for high resolution")

  # Transform scHi-C data (e.g.300kb)
  tflamingo.linear_transformation(code_path = code_path, chr_name=chr_name,resolution=low_res,input_path=paste0(outputs_folder,'/lowres_contact_maps'),  opt_path=paste0(outputs_folder,'/lowres_contact_maps_transformed'),assembly=assembly)
  print("Finished linear transformation for low resolution")
  # Transform scHi-C data (e.g.30kb)
  tflamingo.linear_transformation(code_path = code_path, chr_name=chr_name,resolution=high_res,input_path=paste0(outputs_folder,'/highres_contact_maps'),  opt_path=paste0(outputs_folder,'/highres_contact_maps_transformed'),assembly=assembly)
  print("Finished linear transformation for high resolution")
}

