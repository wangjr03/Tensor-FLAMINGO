library(tFLAMINGOr)

args <- commandArgs(T)

data_path = args[1]
chr_name = args[2]
low_resolution = as.integer(args[3])
high_resolution = as.integer(args[4])
assembly = args[5]
output_dir = args[6]

# first generate contact maps at low-resolution 
tflamingo.generate_matrix(input_path=data_path, resolution=low_resolution, opt_path=paste0(output_dir,'/low_resolution_contact_maps'),assembly=assembly)

# then generate contact maps at high-resolution
tflamingo.generate_matrix(input_path=data_path, resolution=high_resolution, opt_path=paste0(output_dir,'/high_resolution_contact_maps'),assembly=assembly)

# Transform scHi-C data at low resolution
tflamingo.linear_transformation(chr_name=chr_name,resolution=low_resolution,input_path=paste0(output_dir,'/low_resolution_contact_maps'),  opt_path=paste0(output_dir,'/low_resolution_contact_maps_transformed'),assembly=assembly)

# Transform scHi-C data at high resolution
tflamingo.linear_transformation(chr_name=chr_name,resolution=high_resolution,input_path=paste0(output_dir,'/high_resolution_contact_maps'),  opt_path=paste0(output_dir,'/high_resolution_contact_maps_transformed'),assembly=assembly)