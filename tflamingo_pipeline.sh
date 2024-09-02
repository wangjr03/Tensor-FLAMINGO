#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: bash $0 --input_folder <INPUT_FOLDER> --chr_name <CHR_NAME> --low_res <LOW_RES> --high_res <HIGH_RES> --assembly <ASSEMBLY> --outputs_folder <OUTPUTS_FOLDER> --code_path <CODE_PATH>"
    exit 1
}
# Parse the arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in       
        --input_folder)
            INPUT_FOLDER="$2"
            shift
            ;;       
        --chr_name)
            CHR_NAME="$2"
            shift
            ;;
        --low_res)
            LOW_RES="$2"
            shift
            ;;
        --high_res)
            HIGH_RES="$2"
            shift
            ;;
        --assembly)
            ASSEMBLY="$2"
            shift
            ;;
        --outputs_folder)
            OUTPUTS_FOLDER="$2"
            shift
            ;;
        --code_path)
            CODE_PATH="$2"
            shift
            ;;
        *)
            echo "Unknown parameter passed: $1"
            usage
            ;;
    esac
    shift
done

# Check if mandatory parameters are set
if [ -z "$INPUT_FOLDER" ] || [ -z "$CHR_NAME" ] || [ -z "$LOW_RES" ] || [ -z "$HIGH_RES" ] || [ -z "$ASSEMBLY" ] || [ -z "$OUTPUTS_FOLDER" ] || [ -z "$CODE_PATH" ]; then
    echo "Error: --input_folder, --chr_name, --low_res, --high_res, --assembly, --outputs_folder, and --code_path are mandatory."
    usage
fi

# Print the parsed values
# echo "Directory: $dir"
# echo "Data Path: $data_path"
# echo "Resolution: $resolution"
# echo "ASD: $asd"

####
echo "Data preprocessing ..."
mkdir -r $OUTPUTS_FOLDER

Rscript ${CODE_PATH}/data_preprocess.R $INPUT_FOLDER $CHR_NAME $LOW_RES $HIGH_RES $ASSEMBLY $OUTPUTS_FOLDER $CODE_PATH

####
echo "Tensor completing ..."
python ${CODE_PATH}/src/Paralized_Low_rank_tensor_completion_FFTW.py -i ${OUTPUTS_FOLDER}/lowres_contact_maps_transformed -o ${OUTPUTS_FOLDER}/LRTC_low_res_contact_maps -s low_resolution -max_iter 150 -n_core 10
python ${CODE_PATH}/src/Extract_matrix_from_LRTC.py -i ${OUTPUTS_FOLDER}/LRTC_low_res_contact_maps/low_resolution.npy -o ${OUTPUTS_FOLDER}/low_res_contact_maps_FLAMINGO
python ${CODE_PATH}/src/Paralized_Low_rank_tensor_completion_FFTW.py -i ${OUTPUTS_FOLDER}/highres_contact_maps_transformed -o ${OUTPUTS_FOLDER}/LRTC_high_res_contact_maps -s high_resolution -max_iter 150 -n_core 10
python ${CODE_PATH}/src/Extract_matrix_from_LRTC.py -i ${OUTPUTS_FOLDER}/LRTC_high_res_contact_maps/high_resolution.npy -o ${OUTPUTS_FOLDER}/high_res_contact_maps_FLAMINGO

####
echo "3D chromatin structures reconstructing ..."

Rscript ${CODE_PATH}/FLAMINGO_reconstruct.R $OUTPUTS_FOLDER $CHR_NAME $LOW_RES $HIGH_RES

echo "Done!"