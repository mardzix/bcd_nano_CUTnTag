# Multimodal chromatin profiling using nanobody-based single-cell CUT&Tag
### Marek Bartosovic, Goncalo Castelo-Branco
# 
Code repository related to preprint:

https://www.biorxiv.org/content/10.1101/2022.03.08.483459v1

# Step1: prepare environment
### conda environment is provided in env/environment.yaml
```angular2html
conda env create -f env/environment.yaml
```
###Additional package dependencies that need to be installed:


```angular2html
R
install.packages(c('argparse','ggplot2','funr','Signac','scales','Seurat','rmarkdown','mclust','GGally','BiocManager','patchwork','markdown','UpSetR','pheatmap','viridis','purrr','Rmagic','devtools','raster'))
BiocManager::install(c('ensembldb','EnsDb.Mmusculus.v79','GenomeInfoDb', 'GenomicRanges', 'IRanges', 'Rsamtools','BiocGenerics','rtracklayer','limma','slingshot','BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'S4Vectors', 'SingleCellExperiment','SummarizedExperiment', 'batchelor', 'Matrix.utils')))
```
### Have seqtk installed in $PATH
https://github.com/lh3/seqtk

### Install papermill for cli for jupyter notebooks
python3 -m pip install papermill

# Step2: Download the data
use fasterq-dump or alternative to download the fastq files

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198467

# Step3: Clone the github repo with analysis code
```
git clone https://github.com/mardzix/bcd_nano_CUTnTag
```

# Step4: Modify config
Change config/config.yaml and specify path to
1. Absolute path to fastq files
2. Specific path to the tmp folder  (Create any folder, e.g. in home) 
3. Specify conda environment to use 
4. Modify path to cellranger reference


# Step5:  Run the pipeline
Change the cluster profile to the specific cluster profile (e.g. slurm etc.)

see more details here: https://github.com/Snakemake-Profiles/slurm
```
snakemake --snakefile code/workflow/Snakefile_single_modality.smk  --cores 16 --profile htcondor -p                                                                              
```