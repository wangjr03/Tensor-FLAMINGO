# Tensor-FLAMINGO: **Tensor**-based **F**ast **L**ow-r**A**nk **M**atrix completion algorithm for reconstruct**IN**g high-resolution 3D **G**enome **O**rganizations
## Gallery
**The 3D structures of chromosome 21 for 14 single cells in 10-kb resolution** based on Dip-C data (GM12878).
![chr21_Dip-C](./predictions/images/chr21_Dip-C.png)

## Summary
Tensor-FLAMINGO aims to accurately reconstruct the 3D chromatin structures for every single cell from super sparse scHi-C contact maps (missing rate > 99.95%). Remarkably, the tensor-completion-based method borrows information across single cells, while preserving the unique structural variations across single cells.

## Introduction
Tensor-FLAMINGO takes scHi-C data of tens to hundreds of single cells as inputs and reconstructs the single-cell 3D chromatin structures. Tensor-FLAMINGO has two major steps. In the first step, the contact maps of all single cells are modeled as a sparse tensor and then completed using the low-rank tensor completion method. This step gives a dense tensor and imputes the missing values of the original scHi-C data. In the second step, the 3D genome structures are reconstructed for each single cell from the completed chromatin contact map using [FLAMINGO](https://github.com/wangjr03/FLAMINGO/).

## Dependencies
The implementation of the algorithm is based on but not necessarily restricted to `python/3.11.5` and `R/4.3.3`

### R packages
`GenomeInfoDb` 

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomeInfoDb")
```


### python packages
`pyfftw`, `scipy`, `numpy`, `pandas`, `joblib` and `ray`.

```
pip install -r requirements.txt
```

## Installation of Tensor-FLAMINGO
The code for the first step is available in the Github: <br>
```
git clone https://github.com/wangjr03/Tensor-FLAMINGO.git
cd Tensor-FLAMINGO
ls
```
The R package for Tensor-FLAMINGO (*tFlamingorLite*) can be installed through Github using *devtools*:<br>
```
install.packages("devtools")
library(devtools)
install_github('wangjr03/Tensor-FLAMINGO/tFlamingorLite')
```
## Input data
The standard sparse-matrix format scHi-C data is accepted
```
chr1 12345 chr1 13456
```
## Input parameters

_input_folder:_ The folder path to input data, can only contain scHi-C data files.\\
chr_name: The desired chromosome, e.g."chr19"
low_res: The domain-level low resolution used for FLAMINGO reconstruction.
high_res: The bin-level high resolution desired for final results.
assembly: The genome assembly version of input data.
outputs_folder: The folder to store all intermediate and final outputs.
code_path: The path of where tFLAMINGOrLite is, i.e."../Tensor-FLAMINGO"


## One command line to run Tensor FLAMINGO

Usage
```
bash tflamingo_pipeline.sh  --input_folder <INPUT_FOLDER> \
                        --chr_name <CHR_NAME> \
                        --low_res <LOW_RES> \
                        --high_res <HIGH_RES> \
                        --assembly <ASSEMBLY> \
                        --outputs_folder <OUTPUTS_FOLDER> \
                        --code_path <CODE_PATH> 
```

Example 
```
bash tflamingo_pipeline.sh --input_folder "../input_scHiC_data" \
                           --chr_name "chr21" \
                           --low_res 3e5 \
                           --high_res 3e4 \
                           --assembly "hg19" \
                           --outputs_folder "../output" \
                           --code_path "../Tensor-FLAMINGO"

```

## Output data format
For each single cell, a data frame with four columns containing the fragment id (the first column) and the 3D coordinates (the other three columns) will be generated and stored in "OUTPUTS_FOLDER/Tensor-FLAMINGO_results"

## Visualize the 3D genome structure using ParaView
Similar to FLAMINGO, Tensor-FLAMINGO predictions can also be visualized using ParaView. To visualize the 3D genome structure using FLAMINGO, the user need to convert the 3D coordinates into a *.vtk* file. In the **FLAMINGOr** package, a `write.vtk` function is provided for such conversion using the command below:<br>
```
write.vtk(points=res[,2:4],lookup_table=rep(1,dim(res)[1]),name='chr21 5kb 3D structure',opt_path='./chr21_10kb.vtk')
```
*Arguments*:<br>

*points*: 3D coordinates predicted by FLAMINGO in the x,y,z format. <br>

*lookup_table*: The annotation of each point, could be labels or scores, i.e. the compartment PC scores.<br>

*name*: output file name annotated within the file.<br>

*opt_path*: output file path including the file name. <br>

## Step-by-step Instruction 
Here we shown an example of reconstructing the single-cell 3D chromosome structures.
```
tFlamingorLite::tflamingo.data_prepare(code_path, input_folder, low_res, high_res, outputs_folder, chr_name, assembly)
```
Applies tensor-completion method
```
python Tensor-FLAMINGO/src/Paralized_Low_rank_tensor_completion_FFTW.py -i './lowres_contact_maps_transformed' -o './LRTC_low_res_contact_maps' -s low_resolution -max_iter 150 -n_core 10
python Tensor-FLAMINGO/src/Paralized_Low_rank_tensor_completion_FFTW.py -i './highres_contact_maps_transformed' -o './LRTC_high_res_contact_maps' -s high_resolution -max_iter 150 -n_core 10
python Tensor-FLAMINGO/src/Extract_matrix_from_LRTC.py -i './LRTC_low_res_contact_maps/low_resolution.npy' -o './low_res_contact_maps_FLAMINGO'
python Tensor-FLAMINGO/src/Extract_matrix_from_LRTC.py -i './LRTC_high_res_contact_maps/high_resolution.npy' -o './high_res_contact_maps_FLAMINGO'
```
Reconstruct the 3D chromatin structures:
```
library(tFlamingorLite)
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
                                                       chr_name = chr_name, nThread = 4, 
                                                       sample_rate = 0.75, lambda = 10, 
                                                       max_dist = 0.01,
                                                       error_threshold = 1e-3,max_iter = 500)
  
  write.table(res_list[[idx]], file = paste0(outputs_folder,"/Tensor-FLAMINGO_results/Cell_",idx,".txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
}
```
