# Example 
We provide an example dataset containing 50 cells from 10X scRNA 2CLCs mouse data. We will go through the process of using MATES to quantify the TE expression in this sample data. The sample data can be downloaded at [here](https://mcgill-my.sharepoint.com/:u:/g/personal/ruohan_wang4_mail_mcgill_ca/EdwyzdHz1KtLr9G_c503mwsB_tvo0zfH4CCcQJ0XvmM4eQ?e=i4C8ev)

## Package Installation
To install MATES, run the following command in terminal:
```sh
git clone https://github.com/mcgilldinglab/MATES.git

## create new enveriomen
conda create -n mates_env python=3.9
conda activate mates_env
## install required package
conda install -c bioconda samtools
pip install pysam
conda install -c bioconda bedtools
pip install pybedtools
## install MATES
cd MATES
pip3 install -r requirements.txt 
pip install .
pip install velocyto
## add env to jupyter notebook
conda install ipykernel
python -m ipykernel install --user --name=mates_env
```
* [Detailed implementation example](https://github.com/mcgilldinglab/MATES/blob/main/Sample/sample_pipeline.ipynb)