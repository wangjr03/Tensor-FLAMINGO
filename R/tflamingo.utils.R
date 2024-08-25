pred_val <- function(x,lm_list){
    df_list <- split(x,floor(exp(x[,4])/5e6)+1)
    res <- list()
    for(i in names(df_list)){
        i <- as.numeric(i)
        lm_model <- lm_list[[ min(i,length(lm_list)) ]]
        tmp_res <- predict(lm_model,df_list[[as.character(i)]][,-c(1,2)])
        res[[i]] <- data.frame(df_list[[as.character(i)]][,c(1,2)],tmp_res)
    }
    return(res)
}

reformat_matrix <- function(x){
   df[,2] <- as.numeric(as.character(df[,2]))
   df[,4] <- as.numeric(as.character(df[,4]))
   return(df)
}

generate_matrix <- function(df,res,chr_size){
   N <- dim(df)[1]
   chr <- df[1,1]
   loc = match(chr,chr_size[,1])
   M = floor(chr_size[loc,2]/res)+1
   mat <- matrix(0,M,M)
   for(i in 1:N){
      id.1 = floor(df[i,2]/res)+1
      id.2 = floor(df[i,3]/res)+1
      mat[id.1,id.2] = mat[id.2,id.1] = mat[id.2,id.1]+1
   }
   return(mat)
}
