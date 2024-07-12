#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: bash $0 --code_dir <CODE_DIR> --data_path <DATA_PATH> --output_dir <OUTPUT_DIR> --chr_name <CHR_NAME> --assembly <ASSEMBLY> --low_resolution <LOW_RESOLUTION> --high_resolution <HIGH_RESOLUTION>" 
    exit 1
}

# Parse the arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --code_dir)
            CODE_DIR="$2"
            shift
            ;;
        --data_path)
            DATA_PATH="$2"
            shift
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift
            ;;
        --chr_name)
            CHR_NAME="$2"
            shift
            ;;
        --assembly)
            ASSEMBLY="$2"
            shift
            ;;
        --low_resolution)
            LOW_RESOLUTION="$2"
            shift
            ;;
        --high_resolution)
            HIGH_RESOLUTION="$2"
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
if [ -z "$CODE_DIR" ] || [ -z "$DATA_PATH" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$CHR_NAME" ] || [ -z "$ASSEMBLY" ] || [ -z "$LOW_RESOLUTION" ] || [ -z "$HIGH_RESOLUTION" ]; then
    echo "Error: --code_dir, --data_path, --output_dir, --chr_name, --assembly, --low_resolution, and --high_resolution are mandatory."
    usage
fi

# Print the parsed values
echo "Directory: $dir"
echo "Data Path: $data_path"
echo "Resolution: $resolution"
echo "ASD: $asd"

####
echo "Data preprocessing ..."
mkdir -r $OUTPUT_DIR
mkdir -r ${OUTPUT_DIR}/preprocess

Rscript ${CODE_DIR}/src/preprocess_pipeline.r $DATA_PATH $CHR_NAME $LOW_RESOLUTION $HIGH_RESOLUTION $ASSEMBLY ${OUTPUT_DIR}/preprocess

####
echo "Tensor completing ..."

python ${CODE_DIR}/src/Paralized_Low_rank_tensor_completion_FFTW.py -i ${OUTPUT_DIR}/preprocess/low_resolution_contact_maps_transformed -o ${OUTPUT_DIR}/preprocess/LRTC_low_resolution_contact_maps -s low_resolution -max_iter 150 -n_core 10
python ${CODE_DIR}/src/Paralized_Low_rank_tensor_completion_FFTW.py -i ${OUTPUT_DIR}/preprocess/high_resolution_contact_maps_transformed -o ${OUTPUT_DIR}/preprocess/LRTC_high_resolution_contact_maps -s high_resolution -max_iter 150 -n_core 10
python ${CODE_DIR}/src/Extract_matrix_from_LRTC.py -i ${OUTPUT_DIR}/preprocess/LRTC_low_resolution_contact_maps/low_resolution.npy -o ${OUTPUT_DIR}/preprocess/low_resolution_contact_maps_FLAMINGO
python ${CODE_DIR}/src/Extract_matrix_from_LRTC.py -i ${OUTPUT_DIR}/preprocess/LRTC_high_resolution_contact_maps/high_resolution.npy -o ${OUTPUT_DIR}/preprocess/high_resolution_contact_maps_FLAMINGO

####
echo "3D chromatin structures reconstructing ..."

Rscript ${CODE_DIR}/src/reconstruct_pipeline.r $OUTPUT_DIR $CHR_NAME $LOW_RESOLUTION $HIGH_RESOLUTIO

echo "Done!"

