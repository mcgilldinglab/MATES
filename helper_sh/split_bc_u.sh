#!/bin/bash
cat $1.txt | while read line
    do
    echo "Start Processing ${line}" >> ./unique_read/split_barcode.log
    mkdir ./unique_read/${line}
    mkdir ./unique_read/${line}/by_barcode
    python scripts/split_bam_by_barcode.py ./unique_read/${line}_uniqueread.bam \
        STAR_Solo/${line}/${line}_Solo.out/Gene/filtered/barcodes.tsv ./unique_read/${line}/by_barcode/ >> ./multi_read/${line}/${line}_multi_splitting.log
    echo "End Processing ${line}" >> ./unique_read/split_barcode.log
    done

cat $1.txt | while read line
    do
    echo "Start Indexing Unique Bam in ${line}" >> ./unique_read/split_barcode.log
    for file in ./unique_read/${line}/by_barcode/*.bam
        do
        samtools sort ${file} -o ${file}
        samtools index ${file}
        done
    echo "End Indexing Unique Bam in ${line}" >> ./unique_read/split_barcode.log
        done