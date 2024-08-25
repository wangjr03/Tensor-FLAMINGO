tflamingo.load_data <- function(input_PD,input_IF,chr_name){
  pd = read.table(input_PD)
  pd = as.matrix(pd)
  input_if = read.table(input_IF)
  input_if = as.matrix(input_if)
  n = dim(pd)[1]
  
  flamingo_obj = new('flamingo',IF=input_if,PD=pd,n_frag=n,chr_name=chr_name)
  return(flamingo_obj)
}