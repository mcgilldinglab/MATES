#!bin/bash
## extract mapped read names
samtools view sample.bam | cut -f 1 | sort | uniq > mapped_read_names.txt
## seperate unmapped reads from raw fastqs
seqkit grep -v -f mapped_read_names.txt sample_1.fastq.gz > unmapped_R1.fastq
seqkit grep -v -f mapped_read_names.txt sample_2.fastq.gz > unmapped_R2.fastq
## realign unmapped reads with required option
STAR --runThreadN 64 --genomeDir path_to_genomeDir \
--readFilesIn unmapped_R1.fastq unmapped_R2.fastq \
--outFileNamePrefix unmapped_ \
--outFilterMultimapNmax 10 \
--outSAMtype BAM SortedByCoordinate
## seperate unique reads from original bam
samtools view -F 4 -bq 255 sample.bam > sample_unique.bam
## merge realigned bam with unique reads bam
samtools merge realigned.bam sample_unique.bam unmapped_Aligned.sortedByCoord.out.bam
## index realigned bam files
samtools index realigned.bam