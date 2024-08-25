structure_assemble <- function(flamingo_high_res_obj, flamingo_backbone_prediction, flamingo_intra_domain_prediction, inf_dist, max_iter){
  
  ###read all the information
  n = nrow(flamingo_backbone_prediction@coordinates)  #number of domains
  backbone_id <- flamingo_backbone_prediction@id  #domain id
  backbone <- as.matrix(flamingo_backbone_prediction@coordinates)  #backbone coordinates
  n_point_per_domain = flamingo_intra_domain_prediction[[4]]@input_n  #domain size (=100)
  
  ### the mean of sub-diagonal of distance matrix of backbone
  scalar <- as.matrix(dist(backbone))
  scalar <- mean(scalar[which(row(scalar)+1 == col(scalar))])
  
  ### read intra-domain structure
  all_points <- list() # all structures in each domain, with invalid point
  valid_id_list <- list() 
  valid_domain <- c()
  
  for (counter in backbone_id) {  ## loop over all domains
    if(!counter %in% names(flamingo_intra_domain_prediction)) next   ## if id in backbone is not an id in intra prediction, skip
    
    tmp_domain_res = flamingo_intra_domain_prediction[[as.character(counter)]]  ## the intra result of this domain id
    tmp_id <- tmp_domain_res@id ## valid id
    p_t <- as.matrix(tmp_domain_res@coordinates) ##all intra-domain points(coordinates) of one domain
    p_t_valid <- p_t[tmp_id,]  ## valid id points
    tmp_center <- backbone[which(backbone_id == counter),]  ##  backbone coord as center point, (watch out for coding)
    
    old_center <- apply(p_t, 2, mean) ##calculate the column mean of p_t
    radius <- max(apply(p_t, 1, function(x){norm(x-old_center,'2')})) ## calculate the max distance between each point and center as radius
    new_loc <- t(apply(p_t, 1, function(x){tmp_center + (x-old_center)/radius*scalar}))  ##new center
    
    ### save structure
    all_points[[as.character(counter)]] <- as.data.frame(new_loc)
    valid_id_list[[as.character(counter)]] <- tmp_id
    valid_domain = c(valid_domain, counter) ##valid domains, valid means exist in intra-prediction
    
  }
  
  total_struc <- backbone[match(valid_domain, backbone_id),]  ## take the valid backbone, valid means exist in intra_prediction
  
  ## get the id list
  id_list = list()
  for (i in 1:length(valid_domain)) { ## loop over index of valid domain
    tmp_id <- (valid_domain[i]-1)*n_point_per_domain + valid_id_list[[i]]
    id_list[[as.character(valid_domain[i])]] = tmp_id  ##valid id list of whole chr
    
  }
  
  ## start_id_list
  start_id_list = sapply(valid_id_list, function(x)x[1])  ##The first valid id of each valid domain
  end_id_list = sapply(valid_id_list, function(x) tail(x,n=1))  ## The last valid id of each valid domain
  
  ####Prepare input data
  #domain 2 to domain 1 direction, center
  direction <- total_struc[2,]-total_struc[1,]
  direction <- as.vector(direction/norm(direction,'2'))
  
  ### domain1 valid tail to domain 1 center direction
  domain_1_tail_coord = all_points[[1]][end_id_list[1],]
  tail_direction <- domain_1_tail_coord - total_struc[1,]
  tail_direction <- as.vector(as.matrix(tail_direction/norm(tail_direction,'2')))
  
  ###rotation matrix to place the first domain
  r_mat <- rotation_matrix(matrix(tail_direction,ncol = 1), as.matrix(direction,nrow = 1))  #??????
  all_points[[1]] <- rotate(all_points[[1]],total_struc[1,],r_mat) #??????
  
  
  #whole pd 
  pd <- tflamingo_high_res_obj@PD
  error = 10
  d = get_dist_vec(all_points, id_list, pd, inf_dist)  
  t = 1e-4
  
  
  #rotation
  g_p <- matrix(0, length(all_points)-1, 3)  ## all zero matrix of size..
  y_s_p <- matrix(0, length(all_points)-1, 3)
  y_e_p <- matrix(0, length(all_points)-1, 3)
  
  error_p = error_change <- 1  ## assign 1 to two variables
  N <- dim(do.call(rbind, all_points))[1]  ## rbind the rows of all lists in all_points
  p_p <- matrix(0,N,3) 
  iter <- 1
  
  while (error > 1e-4 & error_change > 1e-5) {
    y_s <- get_point(all_points, start_id_list)[2:length(valid_domain),] ## valid start bin positions of domain 2:end
    y_e <- get_point(all_points, end_id_list)[1:length(valid_domain)-1,] ## valid end bin positions of domain 1:(end-1)
    d_tild <- evaluate_dist(y_s,y_e)  ## evaluate the distance between one end point and next start point
    g <- 4*(d_tild-d)*(y_s-y_e)
    
    y = g-g_p
    s = y_s-y_s_p
    t_k = (sum(diag(t(s)%*%s)))/(sum(diag(t(s)%*%y)))
    
    y_s_prim <- y_s-t_k*g
    g_p <- g
    y_s_p <- y_s
    
    #update structure
    for (i in 2:length(all_points)) {
      tmp_points <- all_points[[i]]  ## coord of domain i
      tmp_id <- start_id_list[i]    ## valid start point of domain i
      tmp_start <- tmp_points[tmp_id,] ## coord of valie start point of domaini
      tmp_center <- total_struc[i,]  ## center point of domain i
      old_direction <- as.vector(as.matrix(tmp_start-tmp_center))+1e-3 ##I don't know 
      new_direction <- y_s_prim[i-1,]- tmp_center
      r_mat <- rotation_matrix(old_direction,new_direction)
      all_points[[i]] <- rotate(all_points[[i]],total_struc[i,], r_mat)
    }
    
    #update structure
    if(iter==1){
      y_e_prim <- y_e-t_k*g
      y_e_p <- y_e
      for (i in 1:(length(all_points)-1)) {
        tmp_points <- all_points[[i]]
        tmp_id <- end_id_list[i]
        tmp_end <- tmp_points[tmp_id,]
        tmp_center <- total_struc[i,]
        old_direction <- as.vector(as.matrix(tmp_end - tmp_center))+1e-3
        new_direction <- y_e_prim[i,] - tmp_center
        r_mat <- rotation_matrix(old_direction, new_direction)
        all_points[[i]] <- rotate(all_points[[i]],total_struc[i],r_mat)
      }
    }
    iter <- iter+1
    cur_p <- do.call(rbind,all_points)
    all_points_err <- norm(cur_p-p_p,'F')/norm(cur_p,'F')
    p_p <- cur_p
    error <- norm(d_tild-d,'2')/norm(d,'2')
    error_change <- abs(error-error_p)
    error_p <- error
    
    if(all_points_err<0.0001){
      break
    }
    if(iter>max_iter){
      break
    }
    
  }
  
  all_points_valid = data.frame()
  #all_points_whole = data.frame()
  for (i in 1:length(all_points)) {
    all_points_valid = rbind(all_points_valid, all_points[[i]][valid_id_list[[i]],])
  }
  #for (i in 1:length(all_points)) {
  #  all_points_whole = rbind(all_points_whole, all_points[[i]])
  #}
  
  res <- data.frame(frag_id = unlist(id_list), x = all_points_valid[,1], y = all_points_valid[,2], z = all_points_valid[,3])
  #res_whole <- data.frame(frag_id = unlist(1:4400), x = all_points_whole[,1], y = all_points_whole[,2], z = all_points_whole[,3])
  return(res)
}


rotation_matrix = function(x,y)
{
  
  u=x/sqrt(sum(x^2))
  
  v=y-sum(u*y)*u
  v=v/sqrt(sum(v^2))
  
  cost=sum(x*y)/sqrt(sum(x^2))/sqrt(sum(y^2))
  
  sint=sqrt(max(1-cost^2,0));
  #print(cost)
  return(diag(length(x)) - u %*% t(u) - v %*% t(v) +
    cbind(u,v) %*% matrix(c(cost,-sint,sint,cost), 2) %*% t(cbind(u,v)))
}


#### apply rotation
rotate <- function(x,c,r)
{
  
  tmp_c <- t(apply(x,1,function(y){y-c}))
  z = as.matrix(tmp_c) %*% r
  z = t(apply(z,1,function(y){y+c}))
  return(z)
  
}




get_point <- function(all_points,id_list){
  
  # get start point of each domain, valid
  coord = t(mapply(function(x,y) x[y,],all_points,id_list))
  
  coord = matrix(unlist(coord),ncol=3)
  
  return(coord)
}

evaluate_dist <- function(start_point,end_point){
  
  apply(start_point-end_point,1,function(x){norm(x,'2')})
  
}

get_dist_vec <- function(all_points,id_list,pd,inf_dist){
  n <- length(all_points)   ## number of domains
  start_id <- sapply(id_list[2:(n)],function(y){y[1]})     ## the valid true start id of domain 2:n (absolute start id)
  end_id <- sapply(id_list[1:(n-1)],function(y){tail(y,n=1)})   ## the valid true end id of domain 1:(n-1)    (absolute end id)
  dist_vec <- c()
  # observed distance between the adjacent point in two domain
  for(i in 1:(n-1)){
    
    tmp_dist = pd[start_id[i],end_id[i]]
    
    if(is.na(tmp_dist) | tmp_dist == inf_dist){
      # observation not available
      tmp_dist = ave_dist(start_id[i]-end_id[i],pd)
      
    }else if(tmp_dist == 0){
      tmp_dist = 0.001
    }
    
    dist_vec = c(dist_vec,tmp_dist)
  }
  return(dist_vec)
}


ave_dist <- function(r,pd){
  res = pd[row(pd)-col(pd)==r | col(pd)-row(pd)==r]
  return(mean(res[res<Inf],na.rm=T))
}