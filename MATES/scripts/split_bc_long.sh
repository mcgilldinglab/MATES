#!/bin/bash
if [ -f "$1" ] && [ -f "$2" ]; then
    paste "$1" "$2" | while IFS="$(printf '\t')" read -r line1 line2;
    do
    echo "Start Processing ${line1}" >> ./long_read/split_barcode.log
    mkdir ./long_read/${line1}
    mkdir ./long_read/${line1}/by_barcode
    python MATES/scripts/split_bam_by_bc.py "$2" ./long_read/${line1}.bam \
        ${line2} ./long_read/${line1}/by_barcode/ >> ./long_read/${line1}/${line1}_splitting.log
    echo "End Processing ${line1}" >> ./long_read/split_barcode.log
    done
fi

for line in `cat "$1"`
    do
    echo "Start Indexing Unique Bam in ${line}" >> ./unique_read/split_barcode.log
    for file in ./unique_read/${line}/by_barcode/*.bam;
        do
        samtools sort ${file} -o ${file}
        samtools index ${file}
        done
    echo "End Indexing Unique Bam in ${line}" >> ./unique_read/split_barcode.log
    done
