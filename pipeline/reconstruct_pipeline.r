library(tFLAMINGOr)

args <- commandArgs(T)

output_dir = args[1]
chr_name = args[2]
low_resolution = as.integer(args[3])
high_resolution = as.integer(args[4])

n = length(dir( paste0(output_dir,'/preprocess/high_resolution_contact_maps_FLAMINGO') ))
res_list <- list()
for(idx in 1:n){
  input_PD_high = paste0(output_dir,'/preprocess/high_resolution_contact_maps_FLAMINGO/PD_Cell_',idx,'.txt')
  input_IF_high=paste0(output_dir,'/preprocess/high_resolution_contact_maps_FLAMINGO/IF_Cell_',idx,'.txt')
  input_PD_low=paste0(output_dir,'/preprocess/low_resolution_contact_maps_FLAMINGO/PD_Cell_',idx,'.txt')
  input_IF_low=paste0(output_dir,'/preprocess/low_resolution_contact_maps_FLAMINGO/IF_Cell_',idx,'.txt')
  res_list[[i]] <- tflamingo.main_func(
         input_PD_high=input_PD_high,
         input_IF_high=input_IF_high,
         input_PD_low=input_PD_low,
         input_IF_low=input_IF_low,
         domain_res = low_resolution,
         frag_res = high_resolution,
         chr_name=chr_name,
         downsampling_rates=0.75,
         lambda=10,
         max_dist=0.05,
         nThread=20,
         max_iter=500
  )
}

# save
save(res_list, file=paste0(output_dir,'/result.RData'))
