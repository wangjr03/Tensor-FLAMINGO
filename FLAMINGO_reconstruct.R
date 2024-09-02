library(tFlamingorLite)

args <- commandArgs(T)

outputs_folder = args[1]
chr_name = args[2]
low_res = args[3]
high_res = args[4]

n = length(dir(paste0(outputs_folder,'/high_res_contact_maps_FLAMINGO')))/2
res_list <- list()
if(dir.exists(paste0(outputs_folder,"Tensor-FLAMINGO_results"))){
  print('results dir already exist')
}else{
  dir.create(paste0(outputs_folder,"/Tensor-FLAMINGO_results"))
}
for (idx in 1:n) {
  input_PD_high = paste0(outputs_folder,"/high_res_contact_maps_FLAMINGO/PD_Cell_",idx,".txt")
  input_IF_high = paste0(outputs_folder,"/high_res_contact_maps_FLAMINGO/IF_Cell_",idx,".txt")
  input_PD_low = paste0(outputs_folder,"/low_res_contact_maps_FLAMINGO/PD_Cell_",idx,".txt")
  input_IF_low = paste0(outputs_folder,"/low_res_contact_maps_FLAMINGO/IF_Cell_",idx,".txt")
  
  res_list[[idx]] <- tFlamingorLite::tflamingo.main_func(outputs_folder, idx,input_PD_high, 
                                                       input_IF_high, input_PD_low, input_IF_low, 
                                                       domain_res=low_res, frag_res = high_res, 
                                                       chr_name = chr_name, nThread = 20, 
                                                       sample_rate = 0.75, lambda = 10, 
                                                       max_dist = 0.01,
                                                       error_threshold = 1e-3,max_iter = 500)
  
  write.table(res_list[[idx]], file = paste0(outputs_folder,"/Tensor-FLAMINGO_results/Cell_",idx,".txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
}

#save(res_list, file = paste0(outputs_folder,"/results.RData"))
