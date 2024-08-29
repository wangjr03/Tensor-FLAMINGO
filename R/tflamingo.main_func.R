#' tflamingo.main_func
#'
#' Main function to run FLAMINGO part in Tensor-FLAMINGO
#' @param idx Cell index of the cell currently being calculated
#' @param input_PD_high high-resolution input pairwise distance matrix
#' @param input_IF_high high-resolution input interactive frequency matrix
#' @param input_PD_low low-resolution input pairwise distance matrix
#' @param input_IF_low low-resolution input interactive frequency matrix
#' @param domain_res Size of the domains in bps, e.g. 1e6. Try strawr::readHicBpResolutions() to see available resolutions.
#' @param frag_res Size of the fragment in bps, e.g. 1e4. Try strawr::readHicBpResolutions() to see available resolutions.
#' @param chr_name Name of the chromosome, e.g. chr1. Try strawr::readHicChroms() to see available chromosomes.
#' @param normalization Normalization method in .hic file. Try strawr::readHicNormTypes() to see available methods. Could be 'NONE'.
#' @param nThread Number of thread avalable for the reconstruction. Default = 1.
#' @param sample_rate Fraction of available entries in Hi-C to be used during the reconstruction. Default = 0.75.
#' @param lambda Weights for all sampled entries. Default = 10.
#' @param r Weights for distance between consecutive points. Default = 1.
#' @param max_dist Maximum allowed distance betwee two consecutive points. Default = 0.01
#' @param alpha Convertion factor between interaction frequency and pairwise distance. Default = -0.25.
#' @param inf_dist Maximun allowed distance betwee any two points. Default = 2.
#' @param error_threshold Error thresholds for reconstruction. Default = 1e-3.
#' @param max_iter Maximum iterations. Default = 500.
#' @keywords flamingo_main
#' @return A data.frame containing the FLAMINGO predicted 3D structure.
#' @export


tflamingo.main_func = function(idx,input_PD_high,input_IF_high,input_PD_low,input_IF_low, domain_res,frag_res,chr_name,downsampling_rates,lambda,max_dist,nThread,max_iter){
  temp_folder = paste0('./temp_Cell',idx)
  dir.create(temp_folder)
  dir.create(paste0(temp_folder,'/domain_data'))
  dir.create(paste0(temp_folder,'/genomic_loc'))
  print('Loading datasets')
  tflamingo_high_res_obj = tflamingo.load_data(input_PD_high,input_IF_high, chr_name)
  tflamingo_low_res_obj = tflamingo.load_data(input_PD_low,input_IF_low, chr_name)
  print('Dividing domains...')
  tflamingo.divide_domain(tflamingo_high_res_obj, domain_res, frag_res, temp_folder)
  print('Reconstructing backbones...')
  tflamingo_backbone_prediction = tflamingo_backbone(tflamingo_low_res_obj, sample_rate, lambda, r, max_dist, error_threshold, max_iter,inf_dist)
  print('Reconstructing intra-domain structures...')
  tflamingo_intra_domain_prediction = flamingo_domain(temp_folder, sample_rate, lambda, r, alpha, max_dist, error_threshold, max_iter, inf_dist, nThread)
  print('Assembling structures...')
  res = structure_assemble(tflamingo_high_res_obj, tflamingo_backbone_prediction, tflamingo_intra_domain_prediction,inf_dist,max_iter)
  res$chr = chr_name
  res$start = (res$frag_id-1) * frag_res
  res$end = res$frag_id * frag_res
  res = res[,c('chr','start','end','x','y','z')]
  return(res)
}

