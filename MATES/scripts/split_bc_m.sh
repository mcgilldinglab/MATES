#!/bin/bash

# Check if the provided files exist
if [ -f "$1" ] && [ -f "$2" ]; then
    # Process each line of the input files
    paste "$1" "$2" | while IFS="$(printf '\t')" read -r line1 line2; do
        echo "Start Processing ${line1}" >> ./multi_read/split_barcode.log
        
        # Create necessary directories
        mkdir -p ./multi_read/${line1}/by_barcode
        
        # Get the directory of this script
        SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
        
        # Run the Python script
        python "$SCRIPT_DIR/split_bam_by_bc.py" "$3" ./multi_read/${line1}_multireads.bam \
            ${line2} ./multi_read/${line1}/by_barcode/ "$4" >> ./multi_read/${line1}/${line1}_multi_splitting.log
        
        echo "End Processing ${line1}" >> ./multi_read/split_barcode.log
    done
fi

# Index the BAM files
# for line in $(cat "$1"); do
#     echo "Start Indexing Multi Bam in ${line}" >> ./multi_read/split_barcode.log
    
#     for file in ./multi_read/${line}/by_barcode/*.bam; do
#         samtools sort ${file} -o ${file}
#         samtools index ${file}
#     done
    
#     echo "End Indexing Multi Bam in ${line}" >> ./multi_read/split_barcode.log
# done

for line in $(cat "$1"); do
    echo "Start Indexing Multi Bam in ${line}" >> ./multi_read/split_barcode.log
    export line  # Export `line` so it's accessible to parallel subprocesses
    find ./multi_read/${line}/by_barcode/*.bam | parallel -j "$4" --no-notice "samtools sort {} -o {} && samtools index {}"
    echo "End Indexing Multi Bam in ${line}" >> ./multi_read/split_barcode.log
done
