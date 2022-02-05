library(Seurat)
library(Signac)
library(dplyr)
library(argparse)
library(ggplot2)


parser <- ArgumentParser()

parser$add_argument("--matrix", type="character",
                    help="path to the matrix file")
parser$add_argument("--metadata", type="character",
                    help="path to the metadata")
parser$add_argument("--out", type="character",
                    help="path to the output seurat object Rds")
parser$add_argument("--modality", type="character",default="Unknown",
                    help="modality")


args <- parser$parse_args()

m        <- read.table(file=args$matrix,header = TRUE,stringsAsFactors = FALSE)
metadata <- read.table(file=args$metadata)
rownames(m) <- m[,1]
m           <- m[,-1]
m           <- as.matrix(m)

seurat <- CreateSeuratObject(counts = CreateAssayObject(data = m),
                             meta.data = metadata,
                             assay = 'GA')

VariableFeatures(seurat) <- rownames(seurat)

seurat <- seurat %>% RunSVD()
seurat <- seurat %>% RunUMAP(reduction='lsi',dims=2:30)

Idents(seurat) <- seurat$SubClass
seurat$modality <- args$modality

p1 <- DimPlot(seurat,group.by = 'SubClass',label=TRUE)+ NoLegend()
ggsave(filename = gsub('.Rds','_UMAP.png',args$out),plot = p1,width = 6,height = 6)
saveRDS(object = seurat,file = args$out)
