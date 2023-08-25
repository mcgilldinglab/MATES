#!bin/bash
if [ -f "$1" ] && [ -f "$2" ]; then
    for line in `cat "$1"`
        do
        samtools view -F 4 -bq 255 "$2" > ./unique_read/${line}_uniqueread.bam
        samtools index ./unique_read/${line}_uniqueread.bam
        samtools view -h "$2" | grep -vP "(NH:i:1\b)" | samtools view -Sb > ./multi_read/${line}_multireads.bam
        samtools index ./multi_read/${line}_multireads.bam
        done
else
    echo "File '$1' or '$2' not found."
fi