"""
SpaMFC: Spatial Multi-Feature Clustering for Spatial Transcriptomics Subtype Analysis

Requirements:
- Python 3.10+
- scanpy >= 1.10
- numpy >= 1.20
- pandas >= 2.0
- scikit-learn >= 1.0
- scipy >= 1.7
- matplotlib >= 3.5
- seaborn >= 0.12
- pyyaml >= 6.0
- umap-learn >= 0.5
- numba >= 0.56 (for gene correlation analysis)
- tqdm >= 4.60 (for gene correlation analysis)
- psutil >= 5.8 (for gene correlation analysis)
- statsmodels >= 0.13 (for gene correlation analysis)
- gseapy >= 1.0 (for enrichment analysis)
- infercnvpy >= 0.4 (for CNV analysis)
- openpyxl >= 3.0 (for Excel report generation)
"""

from setuptools import setup, find_packages

packages = find_packages(exclude=["tests", "tests.*"])

setup(
    name="SpaMFC",
    version="2.2.0",
    description="Spatial Multi-Feature Clustering for Spatial Transcriptomics Subtype Analysis",
    author="SpaMFC Team",
    author_email="spamfc@example.com",
    url="https://github.com/spamfc/SpaMFC",
    license="MIT",
    packages=packages,
    install_requires=[
        "scanpy>=1.10",
        "numpy>=1.20",
        "pandas>=2.0",
        "scikit-learn>=1.0",
        "scipy>=1.7",
        "matplotlib>=3.5",
        "seaborn>=0.12",
        "pyyaml>=6.0",
        "umap-learn>=0.5",
        "numba>=0.56",
        "tqdm>=4.60",
        "psutil>=5.8",
        "statsmodels>=0.13",
    ],
    extras_require={
        "enrichment": ["gseapy>=1.0"],
        "cnv": ["infercnvpy>=0.4"],
        "niche": ["scNiche @ git+https://github.com/ZJUFanLab/scNiche.git"],
        "report": ["openpyxl>=3.0"],
        "all": [
            "gseapy>=1.0",
            "infercnvpy>=0.4",
            "scNiche @ git+https://github.com/ZJUFanLab/scNiche.git",
            "openpyxl>=3.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "SpaMFC=SpaMFC.cli:main",
        ],
    },
    python_requires=">=3.10",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)