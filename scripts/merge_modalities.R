#!/usr/bin/env Rscript

library(Seurat)
library(Signac)
library(argparse)
library(UpSetR)
library(purrr)

set.seed(1234)

########### Arguments parser
cat("*** Parsing arguments\n")

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo',
                    help="seurat objects to be merged",nargs="+")
parser$add_argument("-o", "--output", type="character", default='foo',
                    help="output file")
args <- parser$parse_args()
######################## End arguments parser
# args <- readRDS('arguments.Rds')
saveRDS(args,'arguments.Rds')

###################### Load data
get_exp_info_from_path <- function(path){
  path <- strsplit(x = path,split = "/")[[1]]
  result <- list()
  result$bin_size <- strsplit(path[length(path)-1],"_")
  result$bin_size <- unlist(lapply(result$bin_size,function(x){x[length(x)]}))
  result$modality <- unlist(strsplit(path[length(path)-3],"_"))[1]
  result$sample   <- path[length(path)-4]
  return(result)
}

#### Load data #####
cat("*** Loading seurat objects\n")
seurat.ls <- lapply(args$input,readRDS)
####################


samples <- unlist(lapply(seurat.ls,function(x){unique(x$sample)}))
# [1] "bcdCT_MB21_02" "bcdCT_MB21_02" "bcdCT_MB21_03" "bcdCT_MB21_03"


modalities <- unlist(lapply(seurat.ls,function(x){unique(x$modality)}))
# [1] "H3K27ac"  "H3K27me3" "H3K27ac"  "H3K27me3"

print(args$input)
print(samples)
print(modalities)


#######################################
# Find common cells accross the sample
cat("*** Finding common cells accross modalities\n")

cells.ls         <- lapply(seurat.ls,colnames)
cells.common.ls <- list()

upset.plots.ls <- list()

for(sample in unique(samples)){
  sample.index      <- grep(sample,samples)
  cells.common.ls[[sample]]  <- setNames(cells.ls[sample.index],modalities[grep(sample,samples)])
  upset.data                 <- UpSetR::fromList(cells.common.ls[[sample]])
  upset.plots.ls[[sample]]   <- UpSetR::upset(upset.data,order.by = 'freq')
  cells.common.ls[[sample]]  <- do.call('c',cells.common.ls[[sample]])
  cells.common.ls[[sample]]  <- table(cells.common.ls[[sample]])
  cells.common.ls[[sample]]  <- names(cells.common.ls[[sample]][cells.common.ls[[sample]] == length(sample.index)])
}

lapply(upset.plots.ls,print)

#####################################
# Filter for common cells only 
cat("*** Filtering seurat object for common cells\n")
print(seurat.ls)

ncells_before <- lapply(seurat.ls,function(x){dim(x)[2]})

for(sample in unique(samples)){
  sample.index <- grep(sample,samples)
  seurat.ls[sample.index] <- lapply(seurat.ls[sample.index],function(x){
    x[,cells.common.ls[[sample]]]
  })
}
ncells_after <- lapply(seurat.ls,function(x){dim(x)[2]})
#####################################
# Rename feature names
cat(paste0("*** After filtering retained ",ncells_after," out of ",ncells_before," cells \n"))

seurat.modality.ls <- list()
metadata.ls        <- list()
matrix.ls          <- list()

for(modality in unique(modalities)){
  print(modality)
  matrix.ls[[modality]] <- list()
  
  # Merge cells with the same modality
  modalities.to.merge               <- grep(modality,modalities)
  seurat.modality.ls[[modality]]    <- purrr::reduce(seurat.ls[modalities.to.merge],merge)
  seurat.modality.ls[[modality]]    <- AddMetaData(seurat.modality.ls[[modality]],metadata = list("modality"=modality))
  
  # Create merged metadata object
  metadata.ls[[modality]]           <- seurat.modality.ls[[modality]]@meta.data
  colnames(metadata.ls[[modality]]) <- paste0('single.',modality,'.',colnames(metadata.ls[[modality]]))
  
  
  for(assay in names(seurat.modality.ls[[modality]]@assays)){
   cat(paste0("*** Renaming features in assay ", assay,"\n"))

    # Rename the feature names
    newnames <- paste0(modality,'-',seurat.modality.ls[[modality]][[assay]]@counts@Dimnames[[1]])
    
    matrix.ls[[modality]][[assay]]             <- seurat.modality.ls[[modality]][[assay]]@counts
    rownames(matrix.ls[[modality]][[assay]])   <- newnames

    # seurat.modality.ls[[modality]][[assay]]@counts@Dimnames[[1]]       <- newnames
    # seurat.modality.ls[[modality]][[assay]]@data@Dimnames[[1]]         <- newnames
    # rownames(seurat.modality.ls[[modality]][[assay]]@meta.features)    <- newnames
  }
}


# Merge the metadata
metadata.merged     <- purrr::reduce(metadata.ls,cbind)

# Find common assays
common.assays <- lapply(seurat.modality.ls,function(x){names(x@assays)})
common.assays <- do.call('c',common.assays)
common.assays <- table(common.assays)
common.assays <- names(common.assays[common.assays == length(seurat.modality.ls)])

# Merge the count matrices
cat("*** Merging count matrices\n") 

matrix.merged <- list()
for(assay in common.assays){
  matrix.to.merge.ls <- lapply(matrix.ls,function(x){
    x[[assay]]
  })
  matrix.merged[[assay]] <- purrr::reduce(.x = matrix.to.merge.ls,.f = rbind)
}


# Create a merged seurat object
cat("*** Creating new seurat object\n")
seurat.modality.ls[['merged']] <- CreateSeuratObject(counts = matrix.merged[[1]],
                                                     project = 'bcdCT',
                                                     assay = names(matrix.merged)[1],
                                                     meta.data = metadata.merged)
# Add the rest of the assays
cat("*** Adding assays to the seurat object\n")
for(assay in common.assays[2:length(common.assays)]){
  seurat.modality.ls[['merged']][[assay]] <- CreateAssayObject(counts = matrix.merged[[assay]])
}

# Add some metadata
seurat.modality.ls[['merged']] <- AddMetaData(seurat.modality.ls[['merged']],
                                               metadata=list('modality'=paste0(c('merge',unique(modalities)),collapse = '_')))

seurat.modality.ls[['merged']] <- AddMetaData(seurat.modality.ls[['merged']],
                                              metadata=list('sample' = seurat.modality.ls[[1]]$sample))


# Export data
cat("*** Exporting data \n")
dir.create(dirname(args$output),recursive = TRUE)
saveRDS(object = seurat.modality.ls,file = paste0(args$output))

pdf(paste0(dirname(args$output),'/upset_plots.pdf'))
lapply(upset.plots.ls,print)
dev.off()
