# sncRNA workflow
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) 

A small ncRNA workflow in Snakemake.



Installation
------------
To install this workflow, clone the repo:

```
git clone https://github.com/sinanugur/sncRNA-workflows.git
cd sncRNA-workflows

```

Quick start example
-------------------
Create a new directory called data and place your FASTQ files or their symbolic links. You need an active Conda installation with Snakemake. You do not have to install any other requirements. This will trigger a workflow run immediately using 15 threads:

```
snakemake -j 15 --use-conda
```

Databases
---------------------
```

```


Output
------
This workflow will generate `results/` directory. 

This directory contains count tables and sample statistics.

Citiation
---------

Funding
---------
