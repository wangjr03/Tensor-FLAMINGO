library(tFlamingorLite)

args <- commandArgs(T)
input_folder = args[1]
chr_name = args[2]
low_res = as.integer(args[3])
high_res = as.integer(args[4])
assembly= args[5]
outputs_folder = args[6]
code_path = args[7]

tFlamingorLite::tflamingo.data_prepare(code_path, input_folder, low_res, high_res, outputs_folder, chr_name, assembly)