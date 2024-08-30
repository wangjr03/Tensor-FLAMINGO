#' tflamingo.main_func
#'
#' Main function to run FLAMINGO part in Tensor-FLAMINGO
#' @param outputs_folder Folder to store intermediate outputs
#' @param idx Cell index of the cell currently being calculated
#' @param input_PD_high high-resolution input pairwise distance matrix
#' @param input_IF_high high-resolution input interactive frequency matrix
#' @param input_PD_low low-resolution input pairwise distance matrix
#' @param input_IF_low low-resolution input interactive frequency matrix
#' @param domain_res Size of the domains in bps, e.g. 1e6. Try strawr::readHicBpResolutions() to see available resolutions.
#' @param frag_res Size of the fragment in bps, e.g. 1e4. Try strawr::readHicBpResolutions() to see available resolutions.
#' @param chr_name Name of the chromosome, e.g. chr1. Try strawr::readHicChroms() to see available chromosomes.
#' @param nThread Number of thread avalable for the reconstruction. Default = 1.
#' @param sample_rate Fraction of available entries in Hi-C to be used during the reconstruction. Default = 0.75.
#' @param lambda Weights for all sampled entries. Default = 10.
#' @param max_dist Maximum allowed distance betwee two consecutive points. Default = 0.01
#' @param error_threshold Error thresholds for reconstruction. Default = 1e-3.
#' @param max_iter Maximum iterations. Default = 500.
#' @keywords flamingo_main
#' @return A data.frame containing the FLAMINGO predicted 3D structure.
#' @export


tflamingo.main_func = function(outputs_folder,
                               idx,
                               input_PD_high,
                               input_IF_high,
                               input_PD_low,
                               input_IF_low,
                               domain_res,
                               frag_res,
                               chr_name,
                               nThread,
                               sample_rate,
                               lambda,
                               max_dist,
                               error_threshold,
                               max_iter){
  setwd(outputs_folder)
  temp_folder = paste0(outputs_folder,'/temp_Cell',idx)
  dir.create(temp_folder)
  #dir.create(paste0(temp_folder,'/domain_data'))
  #dir.create(paste0(temp_folder,'/genomic_loc'))
  print('Loading datasets')
  tflamingo_high_res_obj = tflamingo.load_data(input_PD_high,input_IF_high, chr_name)
  tflamingo_low_res_obj = tflamingo.load_data(input_PD_low,input_IF_low, chr_name)
  print('Dividing domains...')
  flamingo.divide_domain(temp_folder, tflamingo_high_res_obj, domain_res, frag_res)
  print('Reconstructing backbones...')
  tflamingo_backbone_prediction = flamingo.reconstruct_backbone_structure(tflamingo_low_res_obj, sample_rate, lambda,max_dist,nThread)
  print('Reconstructing intra-domain structures...')
  tflamingo_intra_domain_prediction = flamingo.reconstruct_structure(temp_folder,sample_rate, lambda, max_dist, nThread)
  print('Assembling structures...')
  res = flamingo.assemble_structure(tflamingo_backbone_prediction,tflamingo_high_res_obj,tflamingo_intra_domain_prediction,max_iter)
  return(res)
}

