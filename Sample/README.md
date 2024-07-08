# Example

We provide an example dataset containing 50 cells from 10X scRNA 2CLCs mouse data. This example will guide you through the process of using MATES to quantify TE expression in this sample data. The sample data can be downloaded [here](https://mcgill-my.sharepoint.com/:u:/r/personal/ruohan_wang4_mail_mcgill_ca/Documents/sample_data.tar.gz?csf=1&web=1&e=zQHPAi).

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
pip install .
pip install velocyto

# Add environment to Jupyter Notebook
conda install ipykernel
python -m ipykernel install --user --name=mates_env
```
### For a detailed implementation example, refer to [this Jupyter Notebook.](https://github.com/mcgilldinglab/MATES/blob/main/Sample/sample_pipeline.ipynb)
