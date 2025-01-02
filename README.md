# MATES
A Deep Learning-Based Model for Quantifying Transposable Elements in Single-Cell Sequencing Data ([Nature Communications, 2024](https://www.nature.com/articles/s41467-024-53114-7)).
### Citation
If you use `MATES` in your research, please cite `MATES` publication as follows:
>Wang, R., Zheng, Y., Zhang, Z. et al. MATES: a deep learning-based model for locus-specific quantification of transposable elements in single cell. Nat Commun 15, 8798 (2024). https://doi.org/10.1038/s41467-024-53114-7

## Overview
<img title="Model Overview" alt="Alt text" src="/figures/Model-figure-01.png">

Transposable elements (TEs) are crucial for genetic diversity and gene regulation. Current single-cell quantification methods often align multi-mapping reads to either ‘best-mapped’ or ‘random-mapped’ locations and categorize them at the subfamily levels, overlooking the biological necessity for accurate, locus-specific TE quantification. Moreover, these existing methods are primarily designed for and focused on transcriptomics data, which restricts their adaptability to single-cell data of other modalities. To address these challenges, here we introduce MATES, a deep-learning approach that accurately allocates multi-mapping reads to specific loci of TEs, utilizing context from adjacent read alignments flanking the TE locus. When applied to diverse single-cell omics datasets, MATES shows improved performance over existing methods, enhancing the accuracy of TE quantification and aiding in the identification of marker TEs for identified cell populations. This development facilitates the exploration of single-cell heterogeneity and gene regulation through the lens of TEs, offering an effective transposon quantification tool for the single-cell genomics community.
<!-- MATES is a specialized tool designed for precise quantification of transposable elements (TEs) in various single-cell datasets. The workflow consists of multiple stages to ensure accurate results. In the initial phase, raw reads are mapped to the reference genome, differentiating between unique-mapping and multi-mapping reads associated with TE loci. Unique-mapping reads create coverage vectors (V<sub>u</sub>), while multi-mapping reads remain associated with V<sub>m</sub> vectors, both capturing read distribution around TEs. TEs are then divided into bins, either unique-dominant (U) or multi-dominant (M), based on read proportion. An autoEncoder model is employed to create latent embeddings (Z<sub>m</sub>) capturing local read context and is combined with TE family information (T<sub>k</sub>). In the subsequent stage, the obtained embeddings are used to jointly estimate the multi-mapping ratio (α<sub>i</sub>) via a multilayer perceptron. Training the model involves a global loss (L<sub>1</sub> and L<sub>2</sub>) comprising reconstruction loss and read coverage continuity. Trained to predict multi-mapping ratios, the model counts reads in TE regions, enabling probabilistic TE quantification at the single-cell level. MATES enhances cell clustering and biomarker identification by integrating TE quantification with gene expression methods. -->

<!-- With the burgeoning field of single-cell sequencing data, the potential for in-depth TE quantification and analysis is enormous, opening avenues to gain invaluable insights into the molecular mechanisms underpinning various human diseases. MATES furnishes a powerful tool for accurately quantifying and investigating TEs at specific loci and single-cell level, thereby significantly enriching our understanding of complex biological processes. This opens a new dimension for genomics and cell biology research and holds promise for potential therapeutic breakthroughs. -->

## Relesae Note
* Version 0.1.8: Support TE quantification for data without sufficient multi-mapping TE reads.
* Version 0.1.7: Parallelize preprocessing for 10X-format data.
* Version 0.1.6: Add a simple mode for MATES to quantify TE within 3 lines of code. Add a common errors Q&A. 
* Version 0.1.5: Improve the efficiency of splitting BAM files and counting TEs reads.
<!-- * Version 0.1.4: Enhanced the build_reference.py script and the tutorial to build reference genome for species other than Human and Mouse. -->

MATES is actively under development; please feel free to reach out if you encounter any issue.

## Installation

### Installing MATES
To install MATES, you can run the following command:
```sh
# Clone the MATES repository
git clone https://github.com/mcgilldinglab/MATES.git

# Create a new environment
conda create -n mates_env python=3.9
conda activate mates_env

# Install required packages
conda install -c bioconda samtools -y
conda install -c bioconda bedtools -y

# Install MATES
cd MATES
pip install .

# Add environment to Jupyter Notebook
conda install ipykernel
python -m ipykernel install --user --name=mates_env
```
Installation should take only a few minutes. Verify installation:
```python
import MATES
```

## Usage
### Simple mode
Use the all in one `MATES_pipeline`. Please read our [Examples and APIs](https://github.com/mcgilldinglab/MATES/blob/main/example/simple_mode_sample_pipeline.ipynb) for details. 
```python
from MATES import MATES_pipeline
mates = MATES_pipeline(TE_mode, data_mode, sample_list_file, bam_path_file) #set up parameters
mates.preprocessing() #Preprocessing
mates.run() #train model and quantify both subfamily and locus-level TE expression
```

### Advanced mode
Count coverage vector, Determine U/M region, Generate training and prediction data, Train models, Quantify sub-family level TEs, and Quantify locus_level TEs step by step. Please read our [Examples](https://github.com/mcgilldinglab/MATES/blob/main/example) and [tutorials](https://github.com/mcgilldinglab/MATES/tree/main?tab=readme-ov-file#tutorials). 

### Common Q&A
If you encounter errors when using MATES, please read our [common Q&A](https://github.com/mcgilldinglab/MATES/blob/main/common_errors.md). 

## Tutorials
### Customize the reference genome for the species of interest

We have supported automatic human and mouse TE/Gene reference genome creating using `python build_reference.py --species Human/Mouse`. For Arabidopsis thaliana and Drosophila melanogaster, please visit the shared [folder](https://mcgill-my.sharepoint.com/:f:/g/personal/yumin_zheng_mail_mcgill_ca/EpND9w0Rpf5EpEYGsxoV1I0BUz4svIoqnqWwXcDefeEwFA?e=9p10uW) for GTF file, RepeatMaskers file, and example script to create their TE/Gene reference genome. For other species, please refer to the [tutorial of building TE and Gene reference genome](https://github.com/mcgilldinglab/MATES/blob/main/tutorial/reference_downloading.md).

### Walkthrough Example
From loading data to downstream analysis. Please refer to [Example](https://github.com/mcgilldinglab/MATES/blob/main/example) Section for deatils.

### Example scripts for different type of single cell data
* [MATES pipeline on Smart-seq2 scRNA and 10X scRNA/scATAC/Multi-Omics data](https://github.com/mcgilldinglab/MATES/blob/main/tutorial/pipeline_tutorial.ipynb)
* [MATES pipeline on 10X Visium data](https://github.com/mcgilldinglab/MATES/blob/main/example/10X_visium_example/10X_visium_example.ipynb)
* [MATES pipeline on Long Reads data](https://github.com/mcgilldinglab/MATES/blob/main/tutorial/Tut_LongRead.md)

### 10x scRNA-seq dataset
* [MATES downstream analysis on 10X scRNA data](https://github.com/mcgilldinglab/MATES/blob/main/tutorial/scRNA/10X/scRNA_10X_GeneTE_Analysis.ipynb)
### Smart-seq2 scRNA dataset
* [MATES downstream analysis on Smart-seq2 scRNA data (TE only)](https://github.com/mcgilldinglab/MATES/blob/main/tutorial/scRNA/Smart_seq2/scRNA_SmartSeq_TE_Analysis.ipynb) 
* [MATES downstream analysis on Smart-seq2 scRNA data (Gene+TE)](https://github.com/mcgilldinglab/MATES/blob/main/tutorial/scRNA/Smart_seq2/scRNA_SmartSeq_GeneTE_Analysis.ipynb) 
### 10x scATAC-seq dataset
* [MATES downstream analysis on 10X scATAC data](https://github.com/mcgilldinglab/MATES/blob/main/tutorial/scATAC/scATAC_Peak_TE_analysis.ipynb)

## APIs
The MATES contains **six** modules.
```python
import MATES
from MATES import bam_processor
from MATES import data_processor
from MATES import MATES_model
from MATES import TE_quantifier
from MATES import TE_quantifier_LongRead
from MATES import TE_quantifier_Intronic
```
* **bam_processor**
	The bam_processor module efficiently manages input BAM files by partitioning them into sub-BAM files for individual cells, distinguishing unique mapping from multi mapping reads. It also constructs TE-specific coverage vectors, shedding light on read distributions around TE instances at the single-cell level, enabling accurate TE quantification and comprehensive cellular characterization.
###### For simplicity, in **data_mode**, we use **10X** to represent the format of data in which each BAM file contains reads from multiple cells, and **Smart_seq** to represent the type of data where individual BAM files contain reads from only one cell.


```python 
bam_processor.split_count_10X_data(TE_mode, sample_list_file, bam_path_file, bc_path_file, bc_ind='CB', ref_path = 'Default',num_threads=1)
# Parameters
## TE_mode : <str> exclusive or inclusive, represents whether remove TE instances have overlap with gene (for intronic, refer to below section)
## sample_list_file : <str> path to file conatins sample IDs
## bam_path_file : <str> path to file conatins matching bam file address of sample in sample list
## bc_path_file(optional) : <str> path to file contains matching barcodes list address of sample in sample list
## bc_ind:<str> barcode field indicator in bam files, e.g. CB/CR...
## ref_path(optional): <str> TE reference bed file. Only needed for self generated reference, provide path to reference. By default, exclusive have reference 'TE_nooverlap.bed' and inclusive have reference 'TE_full.bed'.
## num_threads(optional):  <int> The number of process. By default it is 1. Increase the number of threads will reduce the running time, but increase the memory overhead. It should not be larger than the total number of CPUs in your device.
```
```python
bam_processor.split_bam_files(data_mode, threads_num, sample_list_file, bam_path_file,bc_ind = None, bc_path_file=None)
# Parameters
## data_mode : <str> 10X or Smart_seq, we encourage you to use bam_processor.split_count_10X_data() for 10X format data.
## threads_num : <int> Only speed up Smart_seq format data processing. 
## sample_list_file : <str> path to file conatins sample IDs
## bam_path_file : <str> path to file conatins matching bam file address of sample in sample list
## bc_ind:<str> barcode field indicator in bam files, e.g. CB/CR...
## bc_path_file(optional) : <str> path to file contains matching barcodes list address of sample in sample list
```
```python
bam_processor.count_coverage_vec(TE_mode, data_mode, threads_num, sample_list_file, ref_path = "Default", bc_path_file=None)
# Parameters
## TE_mode : <str> exclusive or inclusive, represents whether remove TE instances have overlap with gene (for intronic, refer to below section)
## data_mode : <str> 10X or Smart_seq, we encourage you to use bam_processor.split_count_10X_data() for 10X format data.
## threads_num : <int> Only speed up Smart_seq format data processing. 
## sample_list_file : <str> path to file conatins sample IDs
## ref_path(optional): <str> only needed for self generated reference, provide path to reference. By default, exclusive have reference 'TE_nooverlap.csv' and inclusive have reference 'TE_full.csv'.
## bc_path_file(optional) : <str> only needed for 10X data, path to file contains matching barcodes list address of sample in sample list
```
If you want to perform TE quantification on Long Reads data, you can use **bam_processor.split_bam_files** based on your sequencing plantform. **Instead** of using **bam_processor.count_coverage_vec**, use below function:

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
## Links
Interactive MATES web server: <a>https://mates.cellcycle.org</a>.

## Contact
[Yumin Zheng](mailto:yumin.zheng@mail.mcgill.ca), [Ruohan Wang](mailto:ruohan_wang@brown.edu), [Tao Wu](mailto:tao.wu@bcm.edu), [Jun Ding](mailto:jun.ding@mcgill.ca)
