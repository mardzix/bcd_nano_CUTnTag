#!/usr/bin/env Rscript

library(Seurat)
library(Signac)
library(argparse)
library(purrr)
set.seed(1234)

########### Arguments parser

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo', 
                    help="seurat objects to be merged",nargs="+")
parser$add_argument("-o", "--output", type="character", default='foo', 
                    help="output file")
args <- parser$parse_args()
######################## End arguments parser

# Load the data
seurat.ls <- lapply(args$input,function(x){
  seurat <- readRDS(file=x)
  seurat$orig_file <- x
  return(seurat)
})

# Functions
merge_non_chromatin_assays <- function(assays.ls){
  assay.merged <- purrr::reduce(assays.ls,merge)
  return(assay.merged)
}

merge_chromatin_assays <- function(assays.ls){
  cells_all      <- unlist(lapply(assays.ls,Cells))
  fragments      <- unlist(lapply(assays.ls,function(x){x@fragments}),recursive = TRUE)
  genome         <- unique(unlist(lapply(assays.ls,function(x){unique(genome(x))})))
  
  features.gr    <- unlist(lapply(assays.ls,function(x){x@ranges}),recursive = TRUE)
  features.gr    <- unique(purrr::reduce(features.gr,c))
  features.gr    <- GenomicRanges::reduce(features.gr,min.gapwidth=0)
  
  count.matrix   <- FeatureMatrix(fragments = fragments,
                                  features = features.gr,
                                  cells = cells_all)
  result         <- CreateChromatinAssay(counts = count.matrix,
                                         ranges = features.gr,
                                         fragments = fragments,
                                         genome = genome)
  return(result)
}


# Construct placeholder objects
assays.merged <- list()
seurat.merged <- FALSE

# Get all assay names
all_assays    <- unique(unlist(lapply(seurat.ls,Assays),recursive = TRUE))

# Merge metadata
metadata.merged <- lapply(seurat.ls,function(x){x@meta.data})
metadata.merged <- purrr::reduce(metadata.merged,rbind)

# Merge assay by assay
for(a in all_assays){
  print(a)
  assays.ls    <- lapply(seurat.ls,function(x){x[[a]]})
  is.chromatin <- lapply(assays.ls,function(x){is(x,'ChromatinAssay')})
  print(is.chromatin)
  
  if(sum(unlist(is.chromatin)) < length(is.chromatin)){
    # Not chromatin assay
    assays.merged[[a]] <- merge_non_chromatin_assays(assays.ls)
  } else {
    # Is chromatin assay
    assays.merged[[a]] <- merge_chromatin_assays(assays.ls)
  }
}

# Construct a seurat object
for(a in all_assays){
  if(!is(seurat.merged,'Seurat')){
    seurat.merged <- CreateSeuratObject(counts = assays.merged[[a]],
                                        project = 'bcdCT',
                                        assay = a,meta.data = metadata.merged)
  } else {
    seurat.merged[[a]] <- assays.merged[[a]]    
  }
}

# Export
saveRDS(object = seurat.merged,file = args$output)






