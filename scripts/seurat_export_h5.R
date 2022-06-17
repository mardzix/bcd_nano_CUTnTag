library(Seurat)
library(Signac)
library(argparse)
library(SeuratDisk)



parser <- ArgumentParser()

parser$add_argument("--input", type="character",
                    help="path to the input seurat file")
parser$add_argument("--output", type="character",
                    help="output h5")
args <- parser$parse_args()

# args <- list()
# args$input  <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/multimodal_data/single_modality/ATAC/seurat/peaks/Seurat_object_clustered_renamed.Rds'
# args$output <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/multimodal_data/single_modality/ATAC/seurat/peaks/h5_export/'

seurat <- readRDS(args$input)
SaveH5Seurat(object=seurat, filename = args$output,overwrite=TRUE)