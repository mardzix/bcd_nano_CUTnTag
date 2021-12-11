library(Seurat)
library(Signac)
library(argparse)


set.seed(1234)




########### Arguments parser
cat("*** Parsing arguments\n")

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo',
                    help="seurat object")
parser$add_argument("-o", "--output_file", type="character", default='foo',
                    help="output folder")
args <- parser$parse_args()

###################

cat("*** Loading seurat object \n")
seurat <- readRDS(args$input)
# seurat <- readRDS(file='/data/proj/GCB_MB/bcd_CT/single-cell/results/single_modality/ATAC/seurat_5000/Seurat_object_clustered_renamed.Rds')

cat("*** Finding markers \n")
markers <- FindAllMarkers(seurat)

markers          <- markers[markers$p_val_adj < 0.05,]
markers.positive <- markers[markers$avg_log2FC > 0,  ]

cat("*** Export markers as csv \n")
write.csv(x = markers,file = args$output_file)
write.csv(x = markers,file = gsub(pattern = '\\.csv','_positive.csv',args$output_file))