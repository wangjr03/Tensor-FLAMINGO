#' tflamingo.generate_matrix
#'
#' Generate the contact matrix from scHi-C contact list
#' @param input_path Path containing the scHi-C IF matrices of all chromosomes; will search for the IF matrices in the sub-folder named by the chr_name
#' @param resolution Resolution of the IF matrices, e.g. 10000.
#' @param opt_path path to write the imputed IF matrices.
#' @param assembly Assembly of the scHi-C experiment, e.g. hg38, hg19 and mm10.
#' @keywords tFLAMINGO
#' @return None
#' @export

tflamingo.generate_matrix <- function(input_path,resolution,opt_path,assembly){
    setwd(input_path)
    file_list=dir(input_path)
    suffix=1
    for(fn in file_list){
        tflamingo.process_one_file(fn,resolution,opt_path,assembly,suffix)
        suffix = suffix+1
    }
}