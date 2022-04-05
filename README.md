# sncRNA workflow
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) 

Introduction
------------

This small RNA sequencing pipeline provides a bioinformatics solution to process RNA sequencing data for downstream analyses. The pipeline is built in Snakemake and can be run on different platforms and high performance computing (HPC) systems. It is packed with an Anaconda installation to ensure reproducibility. 



Installation
------------
To install this workflow, clone the repo:

```
git clone https://github.com/sinanugur/sncRNA-workflow.git
cd sncRNA-workflow

```

If you have Anaconda, a new environment can be created

```
conda env create --file environment.yml
conda activate smrnaworkflow

```

Make sure you have human genome file into the `databases` folder, you can download it by typing:

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```

then put this file under `databases` directory. Do not unzip the file, this will also save space.

Quick start example
-------------------
You should create a new directory called __data__ and place your FASTQ files or their symbolic links into `data/` directory. You need an active Conda installation with Snakemake or create an environment using the YAML file. You do not have to install any other requirements. This will trigger a workflow run immediately using 15 threads:

```
snakemake -j 15 --use-conda
```

or if the environment is properly set with all the packages, simply type:

```
snakemake -j 15 
```


Databases
---------------------
```
GENCODE (v38)
miRBase (v22)
piRBase (v1.0)
```
You may update the databases. These versions are the release versions.

Output
------
This workflow will generate `results/` directory. 

This directory contains count tables and sample statistics.

```
Full GENCODE count table: gencode.tsv
miRBase miRNA: miRNA.tsv
miRBase miRNA precursor: miRNA_precursor.tsv
long non-coding RNA (GENCODE): lincRNA.tsv
miscellaneous RNA (GENCODE): misc_RNA.tsv
piwi-interacting RNAs (piRBase): piRNA.tsv
mRNA (GENCODE): protein_coding.tsv
Small nucleolar RNA (GENCODE): snoRNA.tsv
Small nuclear RNA (GENCODE): snRNA.tsv
Small Cajal body-specific RNA (GENCODE): scaRNA.tsv
tRNA (GENCODE): tRNA.tsv
```
Citiation
---------
This workflow was adapted from our small RNA analysis study. Please cite if you find this useful: https://doi.org/10.1080/15476286.2017.1403003

Funding
---------
The study was funded by the European Union’s Horizon 2020 research and innovation program (grant 825741) and the Research Council of Norway under the Program Human Biobanks and Health Data (grant numbers 229621/H10 and 248791/H10). 


