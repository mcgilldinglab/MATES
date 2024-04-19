from setuptools import setup, find_packages
setup(
    name="MATES",
    version="0.1.2",
    packages=find_packages(),
    install_requires=["matplotlib==3.7.2",
        "numpy==1.25.2",
        "pandas==2.0.0",
        "pybedtools==0.9.1",
        "pysam==0.21.0",
        "scipy==1.11.2",
        "torch==2.0.1",
        "tqdm==4.66.1",
        "pyranges==0.0.129"],
)