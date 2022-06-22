#' tflamingo.linear_transformation
#'
#' Transform the all scHi-C IF matrices to the same range as bulk IF with 1D distance-awared linear function
#' @param chr_name Name of the chromosome, e.g. chr1.
#' @param resolution Resolution of the IF matrices, e.g. 10000.
#' @param input_path Path containing the scHi-C IF matrices of all chromosomes; will search for the IF matrices in the sub-folder named by the chr_name
#' @param opt_path path to write the imputed IF matrices.
#' @param assembly Assembly of the scHi-C experiment, e.g. hg38, hg19 and mm10.
#' @keywords tFLAMINGO
#' @return None
#' @export

tflamingo.linear_transformation <- function(chr_name,resolution,input_path,opt_path,assembly){
    chr_size <- GenomicFeatures::getChromInfoFromUCSC(assembly)
    setwd(paste0(input_path,chr_name))
    file <- dir()
    file = file[grep('Cell',file)]
    order_id <- as.numeric(sapply(strsplit(file,'_|[.]'),function(x){x[2]}))
    file <- file[order(order_id,decreasing=F)]
    sc_mat <- list()
    for(i in 1:length(file)){
    sc_mat[[i]] <- as.matrix(as.data.frame(data.table::fread(file[i])))
    }
    # normalize the val towards the averaged count
    averaged_count <- mean(sapply(sc_mat,sum))
    for(i in 1:length(sc_mat)){
        sc_mat[[i]] <- sc_mat[[i]]/sum(sc_mat[[i]])*averaged_count
    }
    # missing rate
    missing_rate <- c()
    for(i in sc_mat){
    missing_rate <- c(missing_rate,1-mean(i>0))
    }
    # given a matrix, get the values, row index and column index
    sc_frame <- list()
    for(i in 1:length(sc_mat)){
        tmp_frame <- Matrix::summary(Matrix::Matrix(sc_mat[[i]],sparse=T))
        id_x <- tmp_frame[,1]
        id_y <- tmp_frame[,2]
        sc_frame[[i]] <- data.frame(x=tmp_frame[,1],y=tmp_frame[,2],dist = abs(tmp_frame[,2]-tmp_frame[,1])*resolution,IF=tmp_frame[,3],
                                    missing_rate=missing_rate[i],if_x = diag(as.matrix(sc_mat[[i]]))[id_x],if_y=diag(as.matrix(sc_mat[[i]]))[id_y])
    }
    pred_sc_if <- list()
    for(i in 1:length(sc_frame)){
        x <- sc_frame[[i]]
        x <- x[,c(1,2,4,3,5)]
        x <- as.matrix(x)
        x <- subset(x,x[,4]>0)
        x[,3] <- log(x[,3])
        x[,4] <- log(x[,4])
        x <- as.data.frame(x)
        x$'ifxy' <- diag(as.matrix(sc_mat[[1]] ))[x[,1]]*diag(as.matrix(sc_mat[[1]] ))[x[,2]]
        transfered_if <- pred_val(x,lm_list)
        transfered_if <- do.call(rbind,transfered_if)
        N <- floor(chr_size[match(chr_name,chr_size[,1]),2]/resolution)+1
        tf_mat <- Matrix::sparseMatrix(
        i = transfered_if[,1], 
        j = transfered_if[,2], 
        x = transfered_if[,3],
        dims = c(N, N), 
        )
        pred_sc_if[[i]] <- as.matrix(exp(tf_mat)-1)
    }
    if(dir.exists(opt_path)){
      print('output dir already exist')
    }else{
      dir.create(opt_path)
    }
    for(i in 1:length(pred_sc_if) ){
        data.table::fwrite(pred_sc_if[[i]],paste0(opt_path,'/RawCount_Cell_',i,".txt"),col.names=F,row.names=F,sep='\t',quote=F,na='NA')
    }

}