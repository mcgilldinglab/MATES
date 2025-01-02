from setuptools import setup, find_packages

setup(
    name="MATES",
    version="0.1.8",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'MATES': ['scripts/*.sh', 'scripts/*.py', 'scripts/Intronic/*.sh', 'scripts/Intronic/*.py']
    },
    install_requires=[
        "matplotlib==3.7.2",
        "numpy==1.25.2",
        "pandas==2.0.0",
        "pybedtools==0.10.0",
        "pysam==0.22.1",
        "scipy==1.11.2",
        "torch==2.0.0",
        "tqdm==4.66.1",
        "pyranges==0.0.129",
        "anndata==0.8.0"
    ],
)
