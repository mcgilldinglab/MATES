# MATES
A Deep Learning-Based Model for Quantifying Transposable Elements in Single-Cell Sequencing Data

## Step 0: Alignment
The raw fastq files are aligned using STAR-Solo for 10X scRNA-seq / scATAC-seq Data and STAR for Smart-Seq2 scRNA-seq Data to reserve multimapping reads. 

- A sample alignment command line for **10X scRNA/scATAC** Data:
```sh
STAR --soloType CB_UMI_Simple --soloCBwhitelist barcode_whitelist \
	--soloMultiMappers EM --runThreadN 64 --genomeDir path_to_genome \
	--outFileNamePrefix STAR_Solo/sample_name --readFilesCommand zcat \
	--readFilesIn sample_name_R1_001.fastq.gz sample_name_R3_001.fastq.gz sample_name_R2_001.fastq.gz \
	--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts \
	--outSAMattributes CR UR CY UY CB UB NH HI
```
> The barcode whitelist can be found at **barcode_whitelist** folder, the original data can be obtained from  https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist- \
> The filtered cell barcodes file provide by STAR-Solo will be at STAR_Solo/sample_name/sample_name.out/Gene/filtered/barcodes.tsv, if you have your own filtered barcodes file, you can simply replace it at the same location.

- A sample alignment command line for **Smaert-seq2 scRNA** Data:
```sh
STAR --runThreadN 64 --genomeDir path_to_genome --readFilesCommand zcat \
        --outFileNamePrefix STAR_Solo/sample/sample_ \
        --readFilesIn sample/sample_1.fastq.gz sample/sample_2.fastq.gz \
        --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
```
## Step 1: Modifying Transposon Element(TE) Reference

The default option of our tool/procedure involves the removal of all transposable element (TE) regions that share overlapping base pairs with gene references. This step is taken to prevent any potential information leakage. To facilitate this process, we provide the TE/Gene references in CSV format. These references can be located within the 'TE_Reference' and 'Gene_Reference' folders. The TE reference data is sourced from https://www.repeatmasker.org/species/mm.html for mouse and https://www.repeatmasker.org/species/hg.html for humans.


