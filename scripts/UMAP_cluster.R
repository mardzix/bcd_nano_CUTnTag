#!/usr/bin/env Rscript

library(Seurat)
library(Signac)
library(ggplot2)
library(argparse)
set.seed(1234)

########### Arguments parser
parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo', 
                    help="seurat object to be clustered")
parser$add_argument("-o", "--output", type="character", default='foo', 
                    help="output file")
parser$add_argument("-a","--assay",type="character",default="bins_5000",
                    help="Assay to be used for dimreduce and clustering")
parser$add_argument('-d','--ndim',type='integer',default=50,
                  help='Number of LSI components to use for UMAP and KNN')
parser$add_argument('-m','--modality',type='character',
                    help='modality to use for clustering if list of objects (only applies for multimodal list of seurat objects)',default='merged')
args <- parser$parse_args()

######################## End arguments parser
# args <- list()
# args$input <- c("results/single_modality/H3K27me3/seurat_5000/Seurat_object.Rds")
# # args$input <- c("results/multiple_modalities/ATAC_H3K27me3/seurat_5000/Seurat_object.Rds")
# args$assay <- 'bins_5000'
# args$ndim <- 50
# args$output <- "/data/proj/GCB_MB/bcd_CT/single-cell/data.Rds"

normalize <- function(seurat_object,assay){
  if (sum(colSums(seurat_object[[assay]]@counts)) == 0){           # Empty counts matrix -> Don't normalize and assume using d
    VariableFeatures(seurat_object) <- rownames(seurat_object)
  }
  else{
    seurat_object <- FindTopFeatures(seurat_object)
    seurat_object <- RunTFIDF(seurat_object)
  }
  return(seurat_object)
}

UMAP_and_cluster <- function(seurat_object, assay, ndim = 50, output = 'seurat_object.Rds'){
  DefaultAssay(seurat_object) <- assay
  if(!'modality' %in% colnames(seurat_object@meta.data)){
    seurat_object$modality <- 'Unknown'
  }
  modality <- unique(seurat_object$modality)
  
  seurat_object <- normalize(seurat_object,assay)
  
  seurat_object <- RunSVD(
    object = seurat_object,
    assay = assay,
    reduction.name = 'lsi'
  )
  
  p.depthcor <- DepthCor(seurat_object)
  ggsave(filename = paste0(dirname(output),'/',modality,'_',assay,'_depthcor.png'),width=4,height=4)
  
  dims          <- c(2:ndim)
  
  seurat_object <- RunUMAP(
    object = seurat_object,
    reduction = 'lsi',
    dims = dims
  )
  
  seurat_object <- FindNeighbors(
    object = seurat_object,
    reduction = 'lsi',
    dims = dims
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
  ggsave(plot = p1+p2,
         filename =  paste0(dirname(output),'/',modality,'_',assay,'_UMAP.png'),width = 8,height = 4)
  return(seurat_object)
 
}

# Load data
seurat.ls <- readRDS(args$input)


# If single modality
if(length(seurat.ls) ==1){
  if (!'modality' %in% colnames(seurat.ls@meta.data)){
    seurat.ls$modality <- 'Unknown'
  }
    seurat.ls <- UMAP_and_cluster(seurat_object = seurat.ls,
                                  assay = args$assay,
                                ndim = args$ndim,
                                output = args$output)
  }

# If multiple modalities
if(length(seurat.ls) > 1){
  seurat.ls[[args$modality]] <- UMAP_and_cluster(seurat_object = seurat.ls[[args$modality]],
                                                 assay = args$assay,
                                                 ndim = args$ndim,
                                                 output = args$output)
}

# Save
saveRDS(object = seurat.ls,file = args$output)