# Example

We provide an example dataset containing 50 cells from 10X scRNA 2CLCs mouse data. This example will guide you through the process of using MATES to quantify TE expression in this sample data. The sample data can be downloaded [here](https://mcgill-my.sharepoint.com/:u:/g/personal/yumin_zheng_mail_mcgill_ca/ESmJrKvPuapBiuh7noTKksIBWosum5rswiG28Sgvb-IDUg?e=5sSfkZ). You may place downloaded data at this folder under the directory with name 'sample_data'.

## Package Installation

To install MATES, run the following commands in the terminal:

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

## Sample workflow template

### For a detailed implementation example, refer to [this Jupyter Notebook.](https://github.com/mcgilldinglab/MATES/blob/main/example/sample_pipeline.ipynb) The run time of this pipeline is about 30 mins. 

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
## build_reference.py requires eight arguments:
## '--species', type=str, choices=['Mouse', 'Human', 'Other'],help='Species type'
## '--ref_mode', type=str, default=['repeats', 'TE'], help='TE reference type'
## '--cut_mode', type=str, default='5prime', choices=['5prime', '3prime'], help='Cut mode'
## '--cut_length', type=int, default=1000, help='Cut length'
## '--intronic', type=bool, default=False, help='Build reference for intronic TE'
## 'other_species_TE', type = str, required=False, help = 'Path to TE reference'
## '--other_species_GTF', type = str, required=False, help = 'Path to GTF of the species'
## '--output_prefix', type = str,required=False, help = 'Output prefix'
python build_reference.py --species Mouse 
python build_reference.py --species Human
```
If you have species other than human/mouse, downloaded TE reference in csv format and gene refrence in GTF format from UCSC table browser, run:
```sh
python build_reference.py --species Other --other_species_TE path_to_TE_reference --other_species_GTF path_to_Gene_reference
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
Please follow the [step by step TE/Gene reference building tutorial](https://github.com/mcgilldinglab/MATES/blob/main/tutorial/reference_downloading.md).

**Example script to generate Drosophila melanogaster TE reference files**

You can follow the the tutorial to download the GTF and RepeatMasker files, or download them from this [shared folder](https://mcgill-my.sharepoint.com/:f:/g/personal/yumin_zheng_mail_mcgill_ca/EirU-9fxxLdCrFYpAbzEOjwBa1TZI9YI4Ery7p3suZoZow?e=KXCHDW). 
```python
python build_reference.py --species Other --other_species_TE Drosophila_TE.csv --other_species_GTF dm6.ensGene.gtf
```
----
### For 10X data
----
### Step 1: Processing Bam Files
To run the first step, you'll be required to furnish three separate .txt files containing essential information: sample names (sample_list_file), respective BAM file addresses (bam_path_file), and in the case of 10X data, the supplementary addresses for barcode files associated with each sample (bc_path_file).
```python
import MATES
from MATES import bam_processor

bam_processor.split_count_10X_data(TE_mode, sample_list_file, bam_path_file, bc_path_file, bc_ind,TE_ref_path, num_threads)
```
#### Alternative solution:
These two commands could take long time to run, but it works well with large bam files.

```python
bam_processor.split_bam_files('10X', threads_num, sample_list_file, bam_path_file, bc_ind, bc_path_file)
bam_processor.count_coverage_vec(TE_mode, '10X', threads_num, sample_list_file, TE_ref_path, bc_path_file)
```

### Step 2: Generate Training/Predicting Samples
```python
from MATES import data_processor

data_processor.calculate_UM_region(TE_mode, '10X', sample_list_file, bin_size=5, proportion=80, ref_path = TE_ref_path, bc_path_file=bc_path_file)

data_processor.generate_training_sample('10X', sample_list_file, bin_size, proportion)

data_processor.generate_prediction_sample(TE_mode, '10X',sample_list_file, bin_size, proportion,ref_path=TE_ref_path, bc_path_file = bc_path_file)
```
### Step 3: Training MATES Model and Make Prediction of α
```python
from MATES import MATES_model

MATES_model.train('10X', sample_list_file, bin_size = 5, proportion = 80, BATCH_SIZE= 4096, AE_LR = 1e-4, MLP_LR = 1e-6, AE_EPOCHS = 200, MLP_EPOCHS = 200, DEVICE = 'cpu')

MATES_model.prediction(TE_mode, '10X', sample_list_file, bin_size = 5, proportion = 80, AE_trained_epochs = 200, MLP_trained_epochs = 200, DEVICE = 'cpu',ref_path = TE_ref_path)
```
### Step 4: Quantify TE Expression Matrix
```python
from MATES import TE_quantifier

TE_quantifier.unique_TE_MTX(TE_mode, '10X', sample_list_file, threads_num, bc_path_file)

TE_quantifier.finalize_TE_MTX('10X', sample_list_file=None)
```

----
For Smart-seq data
----

### Step 1: Processing Bam Files
To run the first step, you'll be required to furnish two separate .txt files containing essential information: sample names, respective BAM file addresses. For Smart-seq data, MATES does not take barcode list files.

```python
import MATES
from MATES import bam_processor

bam_processor.split_bam_files('Smart_seq', threads_num, sample_list_file, bam_path_file, bc_path_file=None)
```
### Step 2: Build Coverage Vectors
```python
bam_processor.count_coverage_vec(TE_mode, 'Smart_seq', threads_num, sample_list_file, bc_path_file=None)
```

### Step 3: Generate Training/Predicting Samples
```python
from MATES import data_processor

data_processor.calculate_UM_region(TE_mode, '10X', sample_list_file, bin_size=5, proportion=80, ref_path = TE_ref_path, bc_path_file=None)

data_processor.generate_training_sample('Smart_seq', sample_list_file, bin_size, proportion)

data_processor.generate_prediction_sample(TE_mode, 'Smart_seq',sample_list_file, bin_size, proportion,ref_path=TE_ref_path, bc_path_file = None)
```
### Step 4: Training MATES Model and Make Prediction of α
```python
from MATES import MATES_model

MATES_model.train('Smart_seq', sample_list_file, bin_size = 5, proportion = 80, BATCH_SIZE= 4096, AE_LR = 1e-4, MLP_LR = 1e-6, AE_EPOCHS = 200, MLP_EPOCHS = 200, DEVICE = 'cpu')

MATES_model.prediction(TE_mode, 'Smart_seq', sample_list_file, bin_size = 5, proportion = 80, AE_trained_epochs = 200, MLP_trained_epochs = 200, DEVICE = 'cpu',ref_path = TE_ref_path)
```
### Step 5: Quantify TE Expression Matrix
```python
from MATES import TE_quantifier

TE_quantifier.unique_TE_MTX(TE_mode, 'Smart_seq', sample_list_file, threads_num, bc_path_file=None)

TE_quantifier.finalize_TE_MTX('Smart_seq', sample_list_file=None)
```
