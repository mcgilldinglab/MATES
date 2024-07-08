<img width="953" alt="image" src="https://github.com/mcgilldinglab/MATES/assets/88182421/a7f1a359-0bc5-46ff-9f1b-f2545a12ac27"># Example

We provide an example dataset containing 50 cells from 10X scRNA 2CLCs mouse data. This example will guide you through the process of using MATES to quantify TE expression in this sample data. The sample data can be downloaded [here](https://mcgill-my.sharepoint.com/:u:/g/personal/ruohan_wang4_mail_mcgill_ca/EdwyzdHz1KtLr9G_c503mwsB6Y4-lawnqEQ1NBx_1Mn9tQ?e=YkcH1B).

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
### For a detailed implementation example, refer to [this Jupyter Notebook.](https://github.com/mcgilldinglab/MATES/blob/main/Sample/sample_pipeline.ipynb)
