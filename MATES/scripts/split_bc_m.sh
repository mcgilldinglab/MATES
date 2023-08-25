#!/bin/bash
if [ -f "$1" ] && [ -f "$2" ]; then
    for line in `cat "$1"`
    do
    echo "Start Processing ${line}" >> ./multi_read/split_barcode.log
    mkdir ./multi_read/${line}
    mkdir ./multi_read/${line}/by_barcode
    python scripts/split_bam_by_bc.py ./multi_read/${line}_multireads.bam \
        "$2" ./multi_read/${line}/by_barcode/ >> ./multi_read/${line}/${line}_multi_splitting.log
    echo "End Processing ${line}" >> ./multi_read/split_barcode.log
    done
fi

for line in `cat "$1"`
    do
    echo "Start Indexing Multi Bam in ${line}" >> ./multi_read/split_barcode.log
    for file in ./multi_read/${line}/by_barcode/*.bam
        do
        samtools sort ${file} -o ${file}
        samtools index ${file}
        done
    echo "End Indexing Multi Bam in ${line}" >> ./multi_read/split_barcode.log
        done