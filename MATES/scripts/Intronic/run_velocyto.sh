#!/bin/bash
#### $1 sample name $2 bamfile $3 barcode file $4 gtfpath

# # Define the conda environment
# CONDA_ENV=your_conda_environment

# # Find the path to conda
# CONDA_EXE=$(which conda)
# CONDA_BASE=$($CONDA_EXE info --base)
# source $CONDA_BASE/etc/profile.d/conda.sh

# # Activate the conda environment
# conda activate $CONDA_ENV

# Check if velocyto is installed, if not install it
if ! command -v velocyto &> /dev/null; then
    echo "velocyto not found, installing..."
    pip install velocyto
else
    echo "velocyto is already installed."
fi

# Create Velocyto directory if it doesn't exist
if [ ! -d Velocyto ]; then
    mkdir -p Velocyto
fi
gtf_path=$(readlink -f "$4")

# Check if all input files exist

if [ -f "$1" ] && [ -f "$2" ] && [ -f "$3" ] && [ -f "$gtf_path" ]; then
    paste "$1" "$2" "$3" | while IFS="$(printf '\t')" read -r line1 line2 line3; do

        # Create sample-specific directory if it doesn't exist
        if [ ! -d "Velocyto/${line1}" ]; then
            mkdir -p "Velocyto/${line1}"
        else
            echo "Directory ${line1} already exists."
        fi

        line2_path=$(readlink -f "${line2}")
        line3_path=$(readlink -f "${line3}")
        echo "Start Running Velocyto for ${line1}"
        cd "Velocyto/${line1}"
        echo "Path to cell barcode file: ${line3_path}"
        echo "Path to bam file: ${line2_path}"
        velocyto run -b "${line3_path}" "${line2_path}" "$gtf_path" --dump p1 > velocyto_progress.log
        echo "End Running Velocyto for ${line1}"
        cd ../..
    done
else
    echo "One or more input files do not exist."
fi

# Deactivate the conda environment
# conda deactivate
