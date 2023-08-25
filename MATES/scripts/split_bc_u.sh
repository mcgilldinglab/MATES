#!/bin/bash
if [ -f "$1" ] && [ -f "$2" ]; then
    for line in `cat "$1"`
    do
    echo "Start Processing ${line}" >> ./unique_read/split_barcode.log
    mkdir ./unique_read/${line}
    mkdir ./unique_read/${line}/by_barcode
    python scripts/split_bam_by_bc.py ./unique_read/${line}_uniqueread.bam \
        "$2" ./unique_read/${line}/by_barcode/ >> ./unique_read/${line}/${line}_unique_splitting.log
    echo "End Processing ${line}" >> ./unique_read/split_barcode.log
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
