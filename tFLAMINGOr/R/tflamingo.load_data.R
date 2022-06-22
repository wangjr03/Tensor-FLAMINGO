#' tflamingo.load_data
#'
#' load the completed scHi-C matrices and make flamingo data obj
#' @param input_PD Path to the completed pair-wise distance matrix.
#' @param input_IF Path to the completed interaction frequency matrix.
#' @param chr_name Name of the chromosome
#' @keywords tFLAMINGO
#' @return A flamingo data object for later use
#' @export

tflamingo.load_data <- function(input_PD,input_IF,chr_name){
    library(FLAMINGOr)
    pd = read.table(input_PD)
    pd = as.matrix(pd)
    input_if = read.table(input_IF)
    input_if = as.matrix(input_if)
    n = dim(pd)[1]
    flamingo_obj = new('flamingo',IF=input_if,PD=pd,n_frag=n,chr_name=chr_name)
    return(flamingo_obj)
}
