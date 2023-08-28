#!/bin/bash
if [ -f "$1" ] && [ -f "$2" ]; then
    paste "$1" "$2" | while IFS="$(printf '\t')" read -r line1 line2;
    do
    echo "Start Processing ${line1}" >> ./multi_read/split_barcode.log
    mkdir ./multi_read/${line1}
    mkdir ./multi_read/${line1}/by_barcode
    python MATES/scripts/split_bam_by_bc.py ./multi_read/${line1}_multireads.bam \
        ${line2} ./multi_read/${line1}/by_barcode/ >> ./multi_read/${line1}/${line1}_multi_splitting.log
    echo "End Processing ${line}" >> ./multi_read/split_barcode.log
    done
fi

for line in `cat "$1"`
    do
    echo "Start Indexing Multi Bam in ${line}" >> ./multi_read/split_barcode.log
    for file in ./multi_read/${line}/by_barcode/*.bam;
        do
        samtools sort ${file} -o ${file}
        samtools index ${file}
        done
    echo "End Indexing Multi Bam in ${line}" >> ./multi_read/split_barcode.log
        done