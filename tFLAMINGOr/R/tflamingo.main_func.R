#' tflamingo.main_func
#'
#' Main function of FLAMINGO, a wraper for all steps
#' @param input_pd_high Path to the high-resolution pair-wise distance matrix.
#' @param input_if_high Path to the high-resolution interaction frequency matrix.
#' @param input_pd_low Path to the low-resolution pair-wise distance matrix.
#' @param input_if_low Path to the low-resolution interaction frequency matrix.
#' @param domain_res Size of the domains in bps, e.g. 300000.
#' @param frag_res Size of the fragment in bps, e.g. 30000.
#' @param chr_name Name of the chromosome, e.g. chr1.
#' @param downsampling_rates Fraction of contacts to be used during the reconstruction.
#' @param lambda Lagrigian coefficient.
#' @param max_dist Maximum allowed distance betwee two consecutive points.
#' @param nThread Number of thread avalable for the reconstruction.
#' @param max_iter Maximum iteration for the assembling algorithm. default 500.
#' @keywords tFLAMINGO
#' @return A data.frame containing the FLAMINGO predicted 3D structure.
#' @export

tflamingo.main_func <- function(input_pd_high,input_if_high,input_pd_low,input_if_low,domain_res,frag_res,chr_name,downsampling_rates,lambda,max_dist,nThread,max_iter=500){
  print('Loading datasets')
  tflamingo_high_res_data_obj = tflamingo.load_data(input_pd_high,input_if_high,chr_name)
  tflamingo_low_res_data_obj = tflamingo.load_data(input_pd_low,input_if_low,chr_name)
  print('Dividing domains...')
  FLAMINGOr::flamingo.divide_domain(flamingo_obj = tflamingo_high_res_data_obj,domain_res=domain_res,frag_res=frag_res)
  print('Reconstructing backbones...')
  tflamingo_backbone_prediction = FLAMINGOr::flamingo.reconstruct_backbone_structure(flamingo_data_obj = tflamingo_low_res_data_obj,sw=downsampling_rates,lambda=lambda,max_dist = max_dist,nThread=1)
  print('Reconstructing intra-domain structures...')
  tflamigo_intra_domain_prediction = FLAMINGOr::flamingo.reconstruct_structure(sw=downsampling_rates,lambda = lambda,max_dist = max_dist,nThread=nThread)
  print('Assembling structures...')
  res = FLAMINGOr::flamingo.assemble_structure(flamingo_backbone_prediction_obj=tflamingo_backbone_prediction,
                              flamingo_final_res_data_obj=tflamingo_high_res_data_obj,
                              list_of_flamingo_domain_prediction_obj=tflamigo_intra_domain_prediction,
                              max_iter=max_iter)
  return(res)
  
}