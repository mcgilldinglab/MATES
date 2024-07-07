# MATES
A Deep Learning-Based Model for Quantifying Transposable Elements in Single-Cell Sequencing Data

## Overview
<img title="Model Overview" alt="Alt text" src="/figures/Model-figure-01.png">
MATES is a specialized tool designed for precise quantification of transposable elements (TEs) in various single-cell datasets. The workflow consists of multiple stages to ensure accurate results. In the initial phase, raw reads are mapped to the reference genome, differentiating between unique-mapping and multi-mapping reads associated with TE loci. Unique-mapping reads create coverage vectors (V<sub>u</sub>), while multi-mapping reads remain associated with V<sub>m</sub> vectors, both capturing read distribution around TEs. TEs are then divided into bins, either unique-dominant (U) or multi-dominant (M), based on read proportion. An autoEncoder model is employed to create latent embeddings (Z<sub>m</sub>) capturing local read context and is combined with TE family information (T<sub>k</sub>). In the subsequent stage, the obtained embeddings are used to jointly estimate the multi-mapping ratio (α<sub>i</sub>) via a multilayer perceptron. Training the model involves a global loss (L<sub>1</sub> and L<sub>2</sub>) comprising reconstruction loss and read coverage continuity. Trained to predict multi-mapping ratios, the model counts reads in TE regions, enabling probabilistic TE quantification at the single-cell level. MATES enhances cell clustering and biomarker identification by integrating TE quantification with gene expression methods.

With the burgeoning field of single-cell sequencing data, the potential for in-depth TE quantification and analysis is enormous, opening avenues to gain invaluable insights into the molecular mechanisms underpinning various human diseases. MATES furnishes a powerful tool for accurately quantifying and investigating TEs at specific loci and single-cell level, thereby significantly enriching our understanding of complex biological processes. This opens a new dimension for genomics and cell biology research and holds promise for potential therapeutic breakthroughs.

## Installation
### Prerequisites
samtools == 1.17
```sh
conda install -c bioconda samtools
pip install pysam
```
bedtools == 2.31.0
```sh
conda install -c bioconda bedtools
pip install pybedtools
```

### Installing MATES
To install MATES, you can run the following command:
```sh
git clone https://github.com/mcgilldinglab/MATES.git
cd MATES
pip3 install -r requirements.txt 
python setup.py install
```
Installation should take only a few minutes. Verify that MATES is correctly installed by running in python:
```sh
import MATES
```
## Links
Interactive MATES web server: <a>https://mates.cellcycle.org</a>.
## Usage
The MATES contains four main modules.
```python
import MATES
from MATES import bam_processor
from MATES import data_processor
from MATES import MATES_model
from MATES import TE_quantifier
```
* **bam_processor**
	The bam_processor module efficiently manages input BAM files by partitioning them into sub-BAM files for individual cells, distinguishing unique mapping from multi mapping reads. It also constructs TE-specific coverage vectors, shedding light on read distributions around TE instances at the single-cell level, enabling accurate TE quantification and comprehensive cellular characterization.
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
```python
bam_processor.count_coverage_vec(TE_mode, data_mode, threads_num, sample_list_file, ref_path = "Default", bc_path_file=None)
# Parameters
## TE_mode : <str> exclusive or inclusive, represents whether remove TE instances have overlap with gene (for intronic, refer to below section)
## data_mode : <str> 10X or Smart_seq
## threads_num : <int>
## sample_list_file : <str> path to file conatins sample IDs
## ref_path(optional): <str> only needed for self generated reference, provide path to reference. By default, exclusive have reference 'TE_nooverlap.csv' and inclusive have reference 'TE_full.csv'.
## bc_path_file(optional) : <str> only needed for 10X data, path to file contains matching barcodes list address of sample in sample list
```
If you want to perform TE quantification on Long Reads data, you can use **bam_processor.split_bam_files** based on your sequencing plantform. **Instead** of using **bam_processor.count_coverage_vec**, use below function:
###### For simplicity, in **data_mode**, we use **10X** to indicating data using barcodes to distinguish data, i.e. you may have a barcode file to seperating the data in the bam file or **Smart_seq** to indicating data do not use barcodes to distinguish data, i.e. one bam file per cell.


```python
bam_processor.count_long_reads(TE_mode, data_mode, threads_num, sample_list_file, bam_dir, ref_path = "Default", bc_path_file=None):
# Parameters
## TE_mode : <str> exclusive or inclusive, represents whether remove TE instances have overlap with gene
## data_mode : <str> 10X or Smart_seq
## threads_num : <int>
## sample_list_file : <str> path to file conatins sample IDs
## bam_dir: <str> path to director conatins sample bam files
## ref_path(optional): <str> only needed for self generated reference, provide path to reference. By default, exclusive have reference 'TE_nooverlap.csv' and inclusive have reference 'TE_full.csv'.
## bc_path_file(optional) : <str> only needed for data using barcodes to distinguish data, path to file contains matching barcodes list address of sample in sample list
```
* **data_processor**
	The data_processor module assists in computing Unique and Multi Regions, generating training samples, and summarizing the expression of multi-mapping reads for prediction.
```python
data_processor.calculate_UM_region(TE_mode, data_mode, sample_list_file, bin_size=5, proportion=80, ref_path = "Default", bc_path_file=None)
# Parameters
## TE_mode : <str> exclusive or inclusive, represents whether remove TE instances have overlap with gene
## data_mode : <str> 10X or Smart_seq
## sample_list_file : <str> path to file conatins sample IDs
## bin_size : <int> size of U/M Region, default = 5
## proportion : <int> proportion of dominated unique reads in U Region / multi reads in M Region, default = 80
## ref_path(optional): <str> only needed for self generated reference, provide path to reference. By default, exclusive have reference 'TE_nooverlap.csv' and inclusive have reference 'TE_full.csv'.
## bc_path_file(optional) : <str> only needed for 10X data, path to file contains matching barcodes list address of sample in sample list
```
```python
data_processor.generate_training_sample(data_mode, sample_list_file, bin_size, proportion)
# Parameters
## data_mode : <str> 10X or Smart_seq
## sample_list_file : <str> path to file conatins sample IDs
## bin_size : <int> size of U/M Region, default = 5
## proportion : <int> proportion of dominated unique reads in U Region / multi reads in M Region, default = 80
```
```python
data_processor.generate_prediction_sample(TE_mode, data_mode,sample_list_file, bin_size, proportion, ref_path = "Default",bc_path_file=None)
# Parameters
## TE_mode : <str> exclusive or inclusive, represents whether remove TE instances have overlap with gene
## data_mode : <str> 10X or Smart_seq
## sample_list_file : <str> path to file conatins sample IDs
## bin_size : <int> size of U/M Region, default = 5
## proportion : <int> proportion of dominated unique reads in U Region / multi reads in M Region, default = 80
## ref_path(optional): <str> only needed for self generated reference, provide path to reference. By default, exclusive have reference 'TE_nooverlap.csv' and inclusive have reference 'TE_full.csv'.
## bc_path_file(optional) : <str> only needed for 10X data, path to file contains matching barcodes list address of sample in sample list
```
* **MATES_model**
	The MATES_model module serves as the core of the MATES framework, encompassing both training and prediction functions. It is responsible for training a neural network model to accurately predict multi-mapping rates of transposable element (TE) instances based on their read coverage vectors.

```python
MATES_model.train(data_mode, sample_list_file, bin_size = 5, proportion = 80, BATCH_SIZE= 4096, AE_LR = 1e-4, MLP_LR = 1e-6, AE_EPOCHS = 200, MLP_EPOCHS = 200, USE_GPU= True)
# Parameters
## data_mode : <str> 10X or Smart_seq
## sample_list_file : <str> path to file conatins sample IDs
## bin_size : <int> size of U/M Region, default = 5
## proportion : <int> proportion of dominated unique reads in U Region / multi reads in M Region, default = 80
## BATCH_SIZE : <int> default = 4096
## AE_LR : <int> learning rate of AutoEncoder, default = 1e-4
## MLP_LR : <int> learning rate of MLP, default = 1e-6
## AE_EPOCHS : <int> training epochs for AutoEncoder, default = 200
## MLP_EPOCHS : <int> training epochs for MLP, default = 200
## USE_GPU : <bool> whether use GU to train the model, default = True
```
```python
MATES_model.prediction(TE_mode, data_mode, sample_list_file, bin_size = 5, proportion = 80, AE_trained_epochs =200, MLP_trained_epochs=200, USE_GPU= True)
# Parameters
## TE_mode : <str> exclusive or inclusive, represents whether remove TE instances have overlap with gene
## data_mode : <str> 10X or Smart_seq
## sample_list_file : <str> path to file conatins sample IDs
## bin_size : <int> size of U/M Region, default = 5
## proportion : <int> proportion of dominated unique reads in U Region / multi reads in M Region, default = 80
## AE_EPOCHS : <int> training epochs for AutoEncoder, default = 200
## MLP_EPOCHS : <int> training epochs for MLP, default = 200
## USE_GPU : <bool> whether use GU to train the model, default = Truet
```
* **TE_quantifier**
	TE_quantifier module facilitates the quantification of TE expression from unique mapping reads and organizes the generation of finalized TE matrix output files.
```python
TE_quantifier.unique_TE_MTX(TE_mode, data_mode, sample_list_file, threads_num, bc_path_file=None)
# Parameters
## TE_mode : <str> exclusive or inclusive, represents whether remove TE instances have overlap with gene
## data_mode : <str> 10X or Smart_seq
## sample_list_file : <str> path to file conatins sample IDs
## threads_num : <int>
## bc_path_file(optional) : <str> only needed for 10X data, path to file contains matching barcodes list address of sample in sample list
```
```python
TE_quantifier.finalize_TE_MTX(data_mode, sample_list_file=None)
# Parameters
## data_mode : <str> 10X or Smart_seq
## sample_list_file(optional) : <str> only needed for 10X data, path to file conatins sample IDs
```
* **TE_quantifier_LongRead**
	TE_quantifier_LongRead module facilitates the quantification of TE expression from unique mapping reads at locus level for Long Read data.
```python
TE_quantifier_LongRead.quantify_locus_TE_MTX(TE_mode, data_mode, sample_list_file)
# Parameters
## TE_mode : <str> exclusive or inclusive, represents whether remove TE instances have overlap with gene
## data_mode : <str> 10X or Smart_seq
## sample_list_file : <str> path to file conatins sample IDs
## long_read : <bool> whether you're quantifying long read data
```
* **TE_quantifier_Intronic**
	TE_quantifier_Intronic module facilitates the quantification of TE expression in Intronic TEs.
```python
implement_velocyto(data_mode, threads_num, sample_list_file, bam_path_file, gtf_path, bc_path_file=None)
# Parameters
## data_mode : <str> 10X or Smart_seq
## threads_num : <int>
## sample_list_file : <str> path to file conatins sample IDs
## bam_path_file : <str> path to file conatins matching bam file address of sample in sample list
## gtf_path : <str> path to the gene gtf file, this is mandatory to implement velocyto
## bc_path_file(optional) : <str> path to file contains matching barcodes list address of sample in sample list
```

```python
parse_velocyto_output(data_mode, threads_num, sample_list_file)
# Parameters
## data_mode : <str> 10X or Smart_seq
## threads_num : <int> threads to use (CPU number)
## sample_list_file : <str> path to file conatins sample IDs
```
```python
count_unspliced_reads(data_mode, threads_num, sample_list_file, ref_path='Default')
# Parameters
## data_mode : <str> 10X or Smart_seq
## threads_num : <int>
## sample_list_file : <str> path to file conatins sample IDs
## ref_path(optional): <str> only needed for self generated reference, provide path to reference. By default TE reference is of name 'TE_intronic.csv'.  
```
```python
count_intornic_coverage_vec(data_mode, threads_num, sample_list_file, ref_path = 'Default',bc_path_file=None)
# Parameters
## data_mode : <str> 10X or Smart_seq
## threads_num : <int> threads to use (CPU number)
## sample_list_file : <str> path to file conatins sample IDs
## ref_path(optional): <str> only needed for self generated reference, provide path to reference. By default TE reference is of name 'TE_intronic.csv'.  
## bc_path_file(optional) : <str> path to file contains matching barcodes list address of sample in sample list
```
```python
generate_prediction_sample(data_mode, sample_list_file, bin_size, proportion, ref_path = 'Default', bc_path_file=None)
# Parameters
## data_mode : <str> 10X or Smart_seq
## sample_list_file : <str> path to file conatins sample IDs
## bin_size : <int> size of U/M Region, default = 5
## proportion : <int> proportion of dominated unique reads in U Region / multi reads in M Region, default = 80 
## ref_path(optional): <str> only needed for self generated reference, provide path to reference. By default TE reference is of name 'TE_intronic.csv'. 
## bc_path_file(optional) : <str> path to file contains matching barcodes list address of sample in sample list
```
```python
quantify_U_TE_MTX(data_mode, sample_list_file)
# Parameters
## data_mode : <str> 10X or Smart_seq
## sample_list_file : <str> path to file conatins sample IDs
```
```python
quantify_M_TE_MTX(data_mode, sample_list_file, bin_size=5, proportion=80, AE_trained_epochs=200, MLP_trained_epochs=200, USE_GPU= True, ref_path = 'Default')
# Parameters
## data_mode : <str> 10X or Smart_seq
## sample_list_file : <str> path to file conatins sample IDs
## bin_size : <int> size of U/M Region, default = 5
## proportion : <int> proportion of dominated unique reads in U Region / multi reads in M Region, default = 80 
## ref_path(optional): <str> only needed for self generated reference, provide path to reference. By default TE reference is of name 'TE_intronic.csv'. 
## AE_EPOCHS : <int> training epochs for AutoEncoder, default = 200
## MLP_EPOCHS : <int> training epochs for MLP, default = 200
## USE_GPU : <bool> whether use GU to train the model, default = Truet
```
```python
correct_intronic_TE(data_mode, sample_list_file, ref_path = 'Default')
# Parameters
## data_mode : <str> 10X or Smart_seq
## sample_list_file : <str> path to file conatins sample IDs
## ref_path(optional): <str> only needed for self generated reference, provide path to reference. By default TE reference is of name 'TE_intronic.csv'. 
``` 
### Step 0: Alignment/TE Reference
#### Alignment
The raw fastq files are aligned using STAR-Solo for 10X scRNA-seq / scATAC-seq Data and STAR for Smart-Seq2 scRNA-seq Data to reserve multimapping reads. 
**Note**: the bam file you input must contains CR and NH fields.

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
> The filtered cell barcodes file provide by STAR-Solo will be at **STAR_Solo/sample_name/sample_name.out/Gene/filtered/barcodes.tsv**, if you have your own filtered barcodes file, you can simply replace it at the same location.

- A sample alignment command line for **Smart-seq2 scRNA** Data:
```sh
STAR --runThreadN 64 --genomeDir path_to_genome --readFilesCommand zcat \
        --outFileNamePrefix STAR_Solo/sample/sample_ \
        --readFilesIn sample/sample_1.fastq.gz sample/sample_2.fastq.gz \
        --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
```
#### TE Reference
We provide two mode of TE reference. **exclusive**, which refers that exclude all TE instances in the reference that have overlapping with gene reference, and **inclusive** refers that we do not remove any TE instances.

The processed TE refernce can be found in TE_ref folder with two different species Human or Mouse. 'TE_nooverlap.csv' is for exclusive mode and 'TE_full.csv' is for inclusive mode. By running the MATES, you will need to place TE refrence file in your working directory.

To generate TE reference, simply run:
```sh
python build_reference.py Mouse
python build_reference.py Human
```
If you have speciese other than human/mouse, downloaded TE reference in csv format and gene refrence in GTF format from UCSC table browser, run:
```sh
python build_reference.py Other path_to_TE_reference path_to_Gene_reference
```
```sh
## A sample of D.melanogaster TE refrence downloaded from UCSC table browser:
$ cat TE_reference.csv | head
#"bin","swScore","milliDiv","milliDel","milliIns","genoName","genoStart","genoEnd","genoLeft","strand","repName","repClass","repFamily","repStart","repEnd","repLeft","id"
"73","845","199","35","16","chr4","130778","131107","-1217024","-","DNAREP1_DM","RC","Helitron","-82","512","268","8"
"74","18658","190","35","41","chr4","1307882","1314210","-33921","-","HETA","LINE","Jockey","-1","6080","1","8"
"585","416","0","0","0","chr4","0","355","-1347776","+","(TTATTATA)n","Simple_repeat","Simple_repeat","1","355","0","8"
"585","15","134","29","29","chr4","688","723","-1347408","+","(TAA)n","Simple_repeat","Simple_repeat","1","35","0","8"
```
You can also find [UCSC table browser use guide](https://github.com/mcgilldinglab/MATES/blob/main/tutorial/reference_downloading.md) to help with downloading reference.

### Step 1: Processing Bam Files
To run the first step, you'll be required to furnish three separate .txt files containing essential information: sample names, respective BAM file addresses, and in the case of 10X data, the supplementary addresses for barcode files associated with each sample.
```python
import MATES
from MATES import bam_processor

bam_processor.split_bam_files(data_mode, threads_num, sample_list_file, bam_path_file, bc_path_file=None)
```
### Step 2: Build Coverage Vectors
```python
bam_processor.count_coverage_vec(TE_mode, data_mode, threads_num, sample_list_file, bc_path_file=None)
```

### Step 3: Generate Training/Predicting Samples
```python
from MATES import data_processor

data_processor.calculate_UM_region(TE_mode, data_mode, sample_list_file, bin_size=5, proportion=80, bc_path_file=None)

data_processor.generate_training_sample(data_mode, sample_list_file, bin_size, proportion)

data_processor.generate_prediction_sample(data_mode,sample_list_file, bin_size, proportion, bc_path_file=None)
```
### Step 4: Training MATES Model and Make Prediction of α
```python
from MATES import MATES_model

MATES_model.train(data_mode, sample_list_file, bin_size = 5, proportion = 80, BATCH_SIZE= 4096, AE_LR = 1e-4, MLP_LR = 1e-6, AE_EPOCHS = 200, MLP_EPOCHS = 200, USE_GPU= True)

MATES_model.prediction(TE_mode, data_mode, sample_list_file, bin_size = 5, proportion = 80, AE_trained_epochs =200, MLP_trained_epochs=200, USE_GPU= True)
```
### Step 5: Quantify TE Expression Matrix
```python
from MATES import TE_quantifier

TE_quantifier.unique_TE_MTX(TE_mode, data_mode, sample_list_file, threads_num, bc_path_file=None)

TE_quantifier.finalize_TE_MTX(data_mode, sample_list_file=None)
```

## Tutorial:
### Pipeline implementation
* [MATES implementation procedures on Smart-seq2 scRNA and 10X scRNA/scATAC/Multi-Omics data](https://github.com/mcgilldinglab/MATES/blob/main/tutorial/pipeline_tutorial.ipynb)
* [MATES implementation procedures on Long Reads data](https://github.com/mcgilldinglab/MATES/blob/main/tutorial/Tut_LongRead.md)

### Training MATES on 10x scRNA-seq dataset
* [MATES downstream analysis on 10X scRNA data](https://github.com/mcgilldinglab/MATES/blob/main/tutorial/scRNA/10X/scRNA_10X_GeneTE_Analysis.ipynb)
### Training MATES on Smart-seq2 scRNA dataset
* [MATES downstream analysis on Smart-seq2 scRNA data (TE only)](https://github.com/mcgilldinglab/MATES/blob/main/tutorial/scRNA/Smart_seq2/scRNA_SmartSeq_TE_Analysis.ipynb) 
* [MATES downstream analysis on Smart-seq2 scRNA data (Gene+TE)](https://github.com/mcgilldinglab/MATES/blob/main/tutorial/scRNA/Smart_seq2/scRNA_SmartSeq_GeneTE_Analysis.ipynb) 
### Training MATES on 10x scATAC-seq dataset
* [MATES downstream analysis on 10X scATAC data](https://github.com/mcgilldinglab/MATES/blob/main/tutorial/scATAC/scATAC_Peak_TE_analysis.ipynb)

