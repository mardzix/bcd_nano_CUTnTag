# Multimodal chromatin profiling using nanobody-based single-cell CUT&Tag
### Marek Bartosovic, Goncalo Castelo-Branco
# 


Code repository related to publication

Bartosovic, M., Castelo-Branco, G. Multimodal chromatin profiling using nanobody-based single-cell CUT&Tag. Nat Biotechnol 41, 794–805 (2023). https://doi.org/10.1038/s41587-022-01535-4


and preprint 


Multimodal chromatin profiling using nanobody-based single-cell CUT&Tag
Marek Bartosovic, Gonçalo Castelo-Branco
bioRxiv 2022.03.08.483459; doi: https://doi.org/10.1101/2022.03.08.483459


# Pipeline 

Standalone pipeline to analyse single-cell nano-CT data is available at:
https://github.com/bartosovic-lab/nanoscope
 

# Data availability
Processed files - Seurat objects, fragments file (cellranger), bigwig tracks per cluster and .h5 matrices are available as supplementary files in the GEO repository
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198467


# Reproducing the analysis
## Step 1: prepare environment
### conda environment is provided in env/environment.yaml
```angular2html
conda env create -f env/environment.yaml
```
### Additional package dependencies that need to be installed:


```angular2html
R
install.packages(c('argparse','ggplot2','funr','Signac','scales','Seurat','rmarkdown','mclust','GGally','BiocManager','patchwork','markdown','UpSetR','pheatmap','viridis','purrr','Rmagic','devtools','raster'))
BiocManager::install(c('ensembldb','EnsDb.Mmusculus.v79','GenomeInfoDb', 'GenomicRanges', 'IRanges', 'Rsamtools','BiocGenerics','rtracklayer','limma','slingshot','BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'S4Vectors', 'SingleCellExperiment','SummarizedExperiment', 'batchelor', 'Matrix.utils')))
```
### Have seqtk installed in $PATH
https://github.com/lh3/seqtk

### Install papermill for cli for jupyter notebooks
```python3 -m pip install papermill```

## Step 2: Download the data
use fasterq-dump or alternative to download the fastq files

```https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump```

GEO repository

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198467


## Step 3: Clone the github repo with analysis code
```
git clone https://github.com/mardzix/bcd_nano_CUTnTag
```

## Step 4: Modify config
Change config/config.yaml and specify path to
1. Absolute path to fastq files
2. Specific path to the tmp folder  (Create any folder, e.g. in home) 
3. Specify conda environment to use 
4. Modify path to cellranger reference


## Step 5:  Run the pipeline
Pipeline is implemented in workflow management software Snakemake 

Change the cluster-specific profile to your preference (e.g. slurm, condor etc.), or run without profile

For some example profiles see: 
- https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
- https://github.com/Snakemake-Profiles/slurm

- https://github.com/Snakemake-Profiles/htcondor
```
snakemake --snakefile code/workflow/Snakefile_single_modality.smk  --cores 16 --profile htcondor -p                                                                              
```

Johannes Köster, Sven Rahmann, Snakemake—a scalable bioinformatics workflow engine, Bioinformatics, Volume 28, Issue 19, 1 October 2012, Pages 2520–2522, https://doi.org/10.1093/bioinformatics/bts480


[![Snakemake](https://img.shields.io/badge/snakemake-≥5.15.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
