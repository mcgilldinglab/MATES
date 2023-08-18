#!/bin/bash
cat $1.txt | while read line
    do
    echo "Start Processing ${line}" >> ./multi_read/split_barcode.log
    mkdir ./multi_read/${line}
    mkdir ./multi_read/${line}/by_barcode
    python scripts/split_bam_by_barcode.py ./multi_read/${line}_multireads.bam \
        STAR_Solo/${line}/${line}_Solo.out/Gene/filtered/barcodes.tsv ./multi_read/${line}/by_barcode/ >> ./multi_read/${line}/${line}_multi_splitting.log
    echo "End Processing ${line}" >> ./multi_read/split_barcode.log
    done

cat $1.txt | while read line
    do
    echo "Start Indexing Multi Bam in ${line}" >> ./multi_read/split_barcode.log
    for file in ./multi_read/${line}/by_barcode/*.bam
        do
        samtools sort ${file} -o ${file}
        samtools index ${file}
        done
    echo "End Indexing Multi Bam in ${line}" >> ./multi_read/split_barcode.log
        done