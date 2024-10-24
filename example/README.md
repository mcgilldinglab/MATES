# Example

We provide an example dataset containing 50 cells from 10X scRNA 2CLCs mouse data. This example will guide you through the process of using MATES to quantify TE expression in this sample data. The sample data can be downloaded [here](https://mcgill-my.sharepoint.com/:u:/g/personal/ruohan_wang4_mail_mcgill_ca/EdwyzdHz1KtLr9G_c503mwsB6Y4-lawnqEQ1NBx_1Mn9tQ?e=YkcH1B). You may place downloaded data at this folder under the directory with name 'sample_data'.

## Package Installation

To install MATES, run the following commands in the terminal:

```sh
# Clone the MATES repository
git clone https://github.com/mcgilldinglab/MATES.git

# Create a new environment
conda create -n mates_env python=3.9
conda activate mates_env

# Install required packages
conda install -c bioconda samtools
pip install pysam
conda install -c bioconda bedtools
pip install pybedtools

# Install MATES
cd MATES
pip3 install -r requirements.txt
python setup.py install
pip install .
pip install velocyto

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
## build_reference.py requires three arguments:
## '--species', type=str, choices=['Mouse', 'Human', 'Other']
## '--cut_mode', type=str, default='5prime', choices=['5prime', '3prime']
## '--cut_length', type=int, default=1000
python build_reference.py --species Mouse 
python build_reference.py --species Human
```
If you have species other than human/mouse, downloaded TE reference in csv format and gene refrence in GTF format from UCSC table browser, run:
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
