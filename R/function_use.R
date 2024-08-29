if2pd <- function(input_if,
                  alpha = -0.25,
                  inf_dist = 2)
{
  
  pd = input_if ^ (alpha)
  pd[which(pd==Inf|is.na(pd))] = inf_dist
  diag(pd) = 0
  
  return(pd)
  
}

pd2gram <- function(pd){
  n = nrow(pd)
  H = Matrix::Diagonal(x=rep(1,n))-1/n*rep(1,n)%*%t(rep(1,n))
  M = -1/2 * H %*% pd^2 %*% H
  M = as.matrix(M)
  
  return(M)
}


get_measurement_set <- function(input_pd, pd_threshold){
  # sparse matrix
  input_pd = as(input_pd,'sparseMatrix')
  col_j = findInterval(seq(input_pd@x)-1,input_pd@p[-1])+1
  row_i = input_pd@i+1
  x_x = input_pd@x
  df = data.frame(row_i,col_j,x_x)
  
  # contact frequency > 0 (if_threshold) ; off-diagonal
  omega = subset(df[,1:2],df[,3] < pd_threshold & df[,1] > df[,2])
  omega = as.matrix(omega)
  
  return(omega)
}

get_bind_set <- function(omega, bind_size){
  
  bind_term = which(omega[,1]-omega[,2] == bind_size)
  n_omega_bind = length(bind_term)
  
  omega_bind = omega[bind_term,,drop=FALSE] # avoid vector with only one sample
  
  return(omega_bind)
  
}

get_sample_set <- function(omega,
                           sample_rate)
{
  n_omega = dim(omega)[1]
  
  # sample 
  sample_set = sample(1:n_omega, sample_rate*n_omega)
  
  omega_sample <<- omega[sample_set,]
  
  return(omega_sample)
  
}

get_element_adjoint_linear <- function(omega){
  
  n_omega = nrow(omega)
  
  # prepare for A*
  w_list = lapply(1:n_omega,function(x){cbind(convert_index(omega[x,]),x)})
  w_list = do.call(rbind,w_list)
  
  all_element = unique(data.frame(w_list[,1:2]))
  loc = prodlim::row.match(data.frame(w_list[,1:2]),all_element)
  
  func_list <- split(data.frame(w_list[,c(3,4)]),loc)
  
  res = list("func_list" = func_list, "all_element" = all_element)
  return(res)
  
}


convert_index <- function(x)
{
  
  df <- expand.grid(as.vector(x),as.vector(x))
  
  df[,3] <- (df[,2] == df[,1])*2-1
  
  return(df)
  
}

linear_proj <- function(omega,x){
  # Linear projection: A(X)
  # fi(X) = <X, omega_i>
  # where i = (i_1, i_2) represents the index of DNA fragment pairs
  # omega_i = e_i1i1 + e_i2i2 - e_i1i2 - e_i2i1
  # where e_ab represents a matrix which has 1 at entry (a,b) and 0 otherwise.
  
  proj <- apply(omega,1,function(y){
    
    y <- as.numeric(y)
    
    p <- x[y[1],y[1]]+x[y[2],y[2]]-x[y[1],y[2]]-x[y[2],y[1]]
    
    return(p)
  })
  
  return(proj)
  
}


flamingo_worker <- function(omega_sample,
                            omega_subdiag,
                            func_list_sample,
                            func_list_subdiag,
                            all_element_sample,
                            all_element_subdiag,
                            b,
                            d,
                            n,
                            lambda,
                            r,
                            error_threshold=1e-3,
                            max_iter=300){
  # initialization
  q = 3
  P_curr <- matrix(rnorm(n*q),n,q)
  P_prev <- matrix(rnorm(n*q),n,q)
  gamma <- rep(0,nrow(omega_sample))
  error <- 10
  
  for (iter in 1:max_iter)
  {
    # calculate gradient
    if(iter == 1){
      grad_prev <- flamingo_grad(P_prev,
                                 omega_sample,
                                 omega_subdiag,
                                 func_list_sample,
                                 func_list_subdiag,
                                 all_element_sample,
                                 all_element_subdiag,
                                 gamma,b,d,lambda,r)
    }else{
      grad_prev <- grad_curr
    }
    
    grad_curr <- flamingo_grad(P_curr,
                               omega_sample,
                               omega_subdiag,
                               func_list_sample,
                               func_list_subdiag,
                               all_element_sample,
                               all_element_subdiag,
                               gamma,b,d,lambda,r)
    
    # BB decent
    P_next = BB_decent(P_prev, P_curr, grad_prev, grad_curr)
    
    # update
    P_prev <- P_curr
    P_curr <- P_next
    error <- norm(P_curr - P_prev,"f")
    
    if (error < error_threshold) break
    
  }
  
  return(P_next)
  
}

flamingo_grad <- function(P,
                          omega_sample,
                          omega_subdiag,
                          func_list_sample,
                          func_list_subdiag,
                          all_element_sample,
                          all_element_subdiag,
                          gamma,b,d,lambda,r)
{
  
  g <- grad_rank(P) +
    lambda * grad_B(P,d,omega_subdiag,func_list_subdiag,all_element_subdiag) +
    r * grad_A(P,gamma,b,omega_sample,func_list_sample,all_element_sample)
  
  return(g)
}

BB_decent <- function(P_prev,
                      P_curr,
                      grad_prev,
                      grad_curr)
{
  
  s = as.matrix(P_curr - P_prev)
  y = as.matrix(grad_curr - grad_prev)
  
  # step size
  t_k = sum(diag( t(s) %*% s )) / sum(diag( t(s) %*% y ))
  
  # update
  P_next = as.matrix(P_curr - t_k * grad_curr)
  
  return(P_next)
  
}

grad_rank <- function(P)
{
  return(2*P)
}

grad_A <- function(P,gamma,b,omega_sample,func_list_sample,all_element_sample){
  
  return(grad_linear(P,-b+gamma,omega_sample,func_list_sample,all_element_sample))
  
}


#### || B(PP^T) - d ||^2
# B^* (B(PP^T) - d) P
grad_B <- function(P,d,omega_subdiag,func_list_subdiag,all_element_subdiag)
{
  
  return(grad_linear(P,-d,omega_subdiag,func_list_subdiag,all_element_subdiag))
  
}


#### || A(PP^T) + y ||^2
# A^* (A(PP^T) + y) P
grad_linear <- function(P,y,omega,func_list,all_element)
{
  
  n = nrow(P)
  x = coord2gram(P) # PP^T
  
  x = linear_proj(omega,x) # A(PP^T)
  
  x = linear_proj_adj(x + y, func_list,all_element,n)
  
  x = x %*% P
  
  return(x)
  
}

coord2gram <- function(P)
{
  
  return(as.matrix(P%*%t(P)))
  
}

linear_proj_adj <- function(x,func_list,all_element,n)
{
  # Linear projection: A(X)
  # fi(X) = <X, omega_i>
  # f*(X) = SUM(X_i * omega_i)
  
  tmp <- sapply(func_list,function(y){
    
    sum(x[y[,2]] * y[,1])
    
  })
  
  mat = Matrix::sparseMatrix(i=all_element[,1],j=all_element[,2],x=tmp,dims=c(n,n))
  
  return(mat)
  
}
