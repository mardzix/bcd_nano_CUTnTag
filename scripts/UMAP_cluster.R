#!/usr/bin/env Rscript

library(Seurat)
library(Signac)
library(ggplot2)
library(argparse)
set.seed(1234)

########### Arguments parser
parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo', 
                    help="seurat objects to be merged",nargs="+")
parser$add_argument("-o", "--output", type="character", default='foo', 
                    help="output file")
parser$add_argument("-a","--assay",type="character",default="bins_5000",
                    help="Assay to be used for dimreduce and clustering")
args <- parser$parse_args()

######################## End arguments parser

# args$input <- c("results/single_modality/H3K27me3/seurat_5000/Seurat_object.Rds")
# args$input <- c("results/multiple_modalities/ATAC_H3K27ac/seurat_5000/Seurat_object.Rds")

seurat.ls <- lapply(args$input,readRDS)

lapply(seurat.ls,function(seurat_object){
  DefaultAssay(seurat_object) <- args$assay
  
  seurat_object <- RunTFIDF(seurat_object)
  seurat_object <- FindTopFeatures(seurat_object)
  
  seurat_object <- RunSVD(
    object = seurat_object,
    assay = args$assay,
    reduction.name = 'lsi'
  )
  
  seurat_object <- RunUMAP(
    object = seurat_object,
    reduction = 'lsi',
    dims = 3:40
  )
  
  seurat_object <- FindNeighbors(
    object = seurat_object,
    reduction = 'lsi',
    dims = 2:40
  )
  
  seurat_object <- FindClusters(
    object = seurat_object,
    algorithm = 3,
    #      resolution = 0.2,
    verbose = FALSE
  )
  
  # seurat_object@meta.data[,paste0('clusters_',antibody)] <- seurat_object@active.ident
  
  p1 <- DimPlot(seurat_object,label=TRUE)
  p2 <- DimPlot(seurat_object,group.by='sample',label=TRUE) + theme(legend.position = 'bottom') + ggtitle(unique(seurat_object$modality))
  p1+p2
  
})