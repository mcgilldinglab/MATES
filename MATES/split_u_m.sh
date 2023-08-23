#!bin/bash
cat $1 | while read line
    do
    samtools view -F 4 -bq 255 ./$2/${line}/${line}_Aligned.sortedByCoord.out.bam > ./unique_read/${line}_uniqueread.bam
    samtools index ./unique_read/${line}_uniqueread.bam
    samtools view -h ./$2/${line}/${line}_Aligned.sortedByCoord.out.bam | grep -vP "(NH:i:1\b)" | samtools view -Sb > ./multi_read/${line}_multireads.bam
    samtools index ./multi_read/${line}_multireads.bam
    done