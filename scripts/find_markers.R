library(Seurat)
library(Signac)


set.seed(1234)




########### Arguments parser
cat("*** Parsing arguments\n")

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo',
                    help="seurat object")
args <- parser$parse_args()

###################

cat("*** Loading seurat object \n")
seurat <- readRDS(args$input)
# seurat <- readRDS(file='/data/proj/GCB_MB/bcd_CT/single-cell/results/single_modality/ATAC/seurat_5000/Seurat_object_clustered_renamed.Rds')

cat("*** Finding markers \n")
markers <- FindAllMarkers(seurat)

cat("*** Export markers as csv \n")
write.csv