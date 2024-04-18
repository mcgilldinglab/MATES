### Step 0: Alignment/TE Reference
#### Alignment
The raw fastq files are aligned using minimap2 or pipeline of your choose.
**Note**: the bam file you input must contains CR and NH fields.

- A sample alignment command line for sample fastq:
```sh
minimap2/minimap2 -ax asm5 index.mmi sample.fastq > long_read_sam/sample.sam
samtools view -bS long_read_sam/sample.sam > long_read_bam/sample.bam
samtools sort long_read_bam/sample.bam -o long_read_bam/sample.bam
samtools index long_read_bam/sample.bam
```
As long reads alignments usually have extremely low multimapping rates, we ignore them at current development.

### Step 1: Split Long Reads by Cell Barcode
```python
bam_processor.split_bam_files(data_mode, threads_num, sample_list_file, bam_path_file,bc_ind = None, bc_path_file=None)
# Parameters
## data_mode : <str> 10X or Smart_seq
## threads_num : <int>
## sample_list_file : <str> path to file conatins sample IDs
## bam_path_file : <str> path to file conatins matching bam file address of sample in sample list
## bc_ind:<str> barcode field indicator in bam files, e.g. CB/CR...
## bc_path_file(optional) : <str> path to file contains matching barcodes list address of sample in sample list
```
### Step 2: Count Locus Specific Reads 
```python
bam_processor.count_long_reads(TE_mode, data_mode, threads_num, sample_list_file, bam_dir, bc_path_file=None)
# Parameters
## TE_mode : <str> exclusive or inclusive, represents whether remove TE instances have overlap with gene
## data_mode : <str> 10X or Smart_seq
## threads_num : <int> threads_num
## sample_list_file : <str> path to file conatins sample IDs
## bam_dir: <str> path to director conatins sample bam files
## bc_path_file(optional) : <str> only needed for data using barcodes to distinguish data, path to file contains matching barcodes list address of sample in sample list
```

### Step 3: Quantify TE Expression Matrix
```python
from MATES import TE_locus_quantifier
unique_locus_TE_MTX(TE_mode, data_mode, sample_list_file, long_read = False, bc_path_file=None)
# Parameters
## TE_mode : <str> exclusive or inclusive, represents whether remove TE instances have overlap with gene
## data_mode : <str> 10X or Smart_seq
## sample_list_file : <str> path to file conatins sample IDs
## bam_dir: <str> path to director conatins sample bam files
## bc_path_file(optional) : <str> only needed for data using barcodes to distinguish data, path to file contains matching barcodes list address of sample in sample list
```
