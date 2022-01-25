cat("*** Loading libraries\n")

library(Seurat)
library(Signac)
library(ggplot2)
library(argparse)

set.seed(1234)

########### Arguments parser
cat("*** Parsing arguments\n")

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo',
                    help="seurat objects to be merged",nargs="+")
parser$add_argument("-o", "--output", type="character", default='foo',
                    help="output file")
parser$add_argument("-m", "--modalities", type="character", default='foo',
                    help="output file",nargs="+")
args <- parser$parse_args()
# 
# args <- list()
# args$input <- c('/data/proj/GCB_MB/bcd_CT/single-cell/results/single_modality/ATAC/seurat_5000/Seurat_object_clustered.Rds',
#                 '/data/proj/GCB_MB/bcd_CT/single-cell/results/single_modality/H3K27ac/seurat_5000/Seurat_object_clustered.Rds',
#                 '/data/proj/GCB_MB/bcd_CT/single-cell/results/single_modality/H3K27me3//seurat_5000/Seurat_object_clustered.Rds')
# args$modalities <- "H3K27ac_H3K27me3"
######## Functions

upset_seurat <- function(seurat.ls,modalities){
  library(UpSetR)
  cat(paste0(c("*** Subseting seurat object with modalities\n",modalities),sep=" "))
  seurat.ls.x <- seurat.ls[modalities]
  
  # Find common samples for the modalities
  samples.common <- lapply(seurat.ls.x,function(x){unique(x$sample)})
  samples.common <- purrr::reduce(samples.common,intersect)
  
  # Keep only cells from common sample
  seurat.ls.x <- lapply(seurat.ls.x,function(x){
    x[,x$sample %in% samples.common]
  })
  
  # List cells 
  cells.ls <- lapply(seurat.ls.x,colnames)
  
  # Upset plot
  upset.data <- UpSetR::fromList(cells.ls)
  p <- UpSetR::upset(upset.data,order.by = 'freq')
  return(p)
}

merge_modalities <- function(seurat.ls,modalities){
  cat(paste0(c("*** Subseting seurat object with modalities",modalities,"\n"),sep=" "))
  seurat.ls.x <- seurat.ls[modalities]
  
  cat(paste0(c("*** Finding common cells accross modalities",modalities,"\n"),sep=" "))
  common.cells <- purrr::reduce(lapply(seurat.ls.x,colnames),intersect)
  
  seurat.ls.x <- lapply(seurat.ls.x,function(x){
    x[,common.cells]
  })
  cat(paste0("*** Finding common assays\n"))
  assays.ls <- lapply(seurat.ls.x,function(x){names(x@assays)})
  common.assays <- purrr::reduce(assays.ls,intersect)
  
  matrix.ls <- list()
  
  # Get list of matrices to merge
  cat(paste0("*** Merge matrices \n"))
  for(assay in common.assays){
    matrix.ls[[assay]] <- lapply(seurat.ls.x,function(x){
      matrix           <- x[[assay]]@counts
      modality         <- unique(x$modality)
      rownames(matrix) <- paste0(modality,'.',rownames(matrix))
      matrix
      })
    matrix.ls[[assay]] <- purrr::reduce(matrix.ls[[assay]],rbind)
  }
  
  # Get and merge metadata
  cat(paste0("*** Merge metadata \n"))
  metadata.ls <- lapply(seurat.ls.x,function(x){
    metadata           <- x@meta.data
    modality           <- unique(x$modality)
    colnames(metadata) <- paste0(modality,'.',colnames(metadata))
    metadata
    })
  
  metadata <- purrr::reduce(metadata.ls,cbind)
  
  cat(paste0("*** Create new seurat object \n"))
  new_obj <- CreateSeuratObject(counts = matrix.ls[[1]],
                                project = 'bcdCT',
                                assay = names(matrix.ls)[1],
                                meta.data = metadata)
  
  cat(paste0("*** Add more modalities \n"))
  for(m in names(matrix.ls[2:length(matrix.ls)])){
    cat(paste0(m,"\n"))
    new_obj[[m]] <- CreateAssayObject(counts = matrix.ls[[m]])
  }
  
  new_obj$modality <- paste0(modalities,collapse="_")
  new_obj$sample   <- seurat.ls.x[[1]]$sample
  
  return(new_obj)
}
##################
cat("*** Loading seurat objects\n")
seurat.ls <- lapply(args$input,readRDS)

seurat.modalities       <- unlist(lapply(seurat.ls, function(x){unique(x$modality)}))
names(seurat.ls)        <- seurat.modalities

modalities.c            <- unlist(strsplit(args$modalities,'_'))

# Create directory
dir.create(dirname(args$output),recursive = TRUE)

# Upset plot of cells overlap
pdf(paste0(dirname(args$output),'/upset_plots.pdf'),width=4,height=4)
upset_seurat(seurat.ls,modalities.c)
dev.off()

# Merge 
seurat.ls[['merged']] <- merge_modalities(seurat.ls,modalities.c)
cat("*** Saving the object \n")
saveRDS(object = seurat.ls,file = args$output)




