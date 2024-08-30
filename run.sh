#!/bin/bash -login
#SBATCH --job-name=Tensor-FLAMINGO          # Job name
#SBATCH --time=04:00:00             # Time limit (hh:mm:ss)
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=4                 # Number of tasks (processes)
#SBATCH --cpus-per-task=4           # Number of CPUs per task
#SBATCH --mem=100G

module purge 
module load R/4.3.3
module load Python/3.11.5

dos2unix /.../Tensor-FLAMINGO/tflamingo_pipeline.sh

bash /.../tflamingo_pipeline.sh --input_folder "/..." --chr_name "chr21" --low_res 3e5 --high_res 3e4 --assembly "hg19" --outputs_folder "/..." --code_path "/.../Tensor-FLAMINGO"



