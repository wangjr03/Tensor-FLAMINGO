flamingo_basic <- function(input_pd, sample_rate, lambda, r, max_dist, error_threshold, max_iter, inf_dist){
  input_pd = as.matrix(input_pd)
  n = nrow(input_pd)
  input_pd[which(input_pd==Inf|is.na(input_pd))] <- inf_dist
  rm_id = which(apply(input_pd,1,min)==inf_dist) #invalid idx
  diag(input_pd)=0
  
  
  M = pd2gram(input_pd)
  ##define measurement set omega
  omega = get_measurement_set(input_pd, inf_dist)
  n_omega = dim(omega)[1]
  
  omega_subdiag = get_bind_set(omega,1)
  n_omega_subdiag = nrow(omega_subdiag)
  
  if (n_omega_subdiag==0) {
    return(NULL)
  }
  
  omega_sample = get_sample_set(omega,sample_rate)
  n_omega_sample = dim(omega_sample)[1]  
  
  
  # prepare for A*
  precal_sample = get_element_adjoint_linear(omega_sample)
  func_list_sample = precal_sample$func_list
  all_element_sample = precal_sample$all_element
  
  # prepare for B*
  precal_subdiag = get_element_adjoint_linear(omega_subdiag)
  func_list_subdiag = precal_subdiag$func_list
  all_element_subdiag = precal_subdiag$all_element
  
  #### pre-calculate the b and d
  b = linear_proj(omega_sample,M)
  
  d = linear_proj(omega_subdiag,M)
  
  # control the sub-diagonal
  for(i in 1:length(d)){
    d[i] = min(d[i],max_dist)
  }
  
  P <- flamingo_worker(omega_sample, 
                       omega_subdiag, 
                       func_list_sample, 
                       func_list_subdiag, 
                       all_element_sample, 
                       all_element_subdiag, 
                       b,d,n,lambda,r,error_threshold,max_iter)
  
  
  #### keep the valid samples
  if(length(rm_id)>0){
    frag_id <- (1:n)[-rm_id]
  }else{
    frag_id <- 1:n
  }
  return(new('flamingo_prediction',id = frag_id, coordinates = P,input_n = n))
}


