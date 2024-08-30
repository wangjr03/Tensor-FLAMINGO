#' tflamingo.process_one_file
#'
#' Generate the contact matrix from scHi-C contact list
#' @param input_path Path containing the scHi-C IF matrices of all chromosomes; will search for the IF matrices in the sub-folder named by the chr_name
#' @param resolution Resolution of the IF matrices, e.g. 10000.
#' @param chr_name chromosome name e.g. "chr19"
#' @param opt_path path to write the imputed IF matrices.
#' @param assembly Assembly of the scHi-C experiment, e.g. hg38, hg19 and mm10.
#' @keywords tFLAMINGO
#' @return None
#' @export

tflamingo.process_one_file <- function(input_path,resolution,chr_name,opt_path,assembly,suffix){
    raw_data <- read.table(input_path)
    raw_data[,1] <- as.character(raw_data[,1])
    raw_data[,3] <- as.character(raw_data[,3])
    raw_data <- subset(raw_data,raw_data[,1] == raw_data[,3])
    chr_size <- getChromInfoFromUCSC(assembly)
    raw_df <- raw_data[,c(1,2,4)]
    chr_split_df <- split(raw_df,raw_df[,1])

    path <- paste0(opt_path,"/",chr_name)
    if(!dir.exists(path)){
      dir.create(path,recursive = TRUE)
    }
    # check if there is any data for ith chr
    if(!chr_name %in% names(chr_split_df)){
        print("no data for this chromosome")
    }
    tmp_mat <- generate_matrix(chr_split_df[[chr_name]],resolution,chr_size)
    write.table(tmp_mat,paste0(path,"/Cell_",suffix,".txt"),col.names=F,row.names=F,sep='\t',quote=F)

}
