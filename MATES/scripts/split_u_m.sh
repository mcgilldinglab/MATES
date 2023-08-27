#!bin/bash
if [ -f "$1" ] && [ -f "$2" ]; then
    paste "$1" "$2" | while IFS="$(printf '\t')" read -r line1 line2;
        do
        samtools view -F 4 -bq 255 "${line2}" > ./unique_read/${line1}_uniqueread.bam
        samtools index ./unique_read/${line1}_uniqueread.bam
        samtools view -h "${line2}" | grep -vP "(NH:i:1\b)" | samtools view -Sb > ./multi_read/${line1}_multireads.bam
        samtools index ./multi_read/${line1}_multireads.bam
        done
else
    echo "File '$1' or '$2' not found."
fi