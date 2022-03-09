library(Seurat)
library(Signac)
library(SeuratWrappers)
library(SeuratDisk)
library(argparse)


set.seed(1234)

########### Arguments parser
cat("*** Parsing arguments\n")

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo',
                    help="seurat object")
parser$add_argument("-c", "--clusters", type="character", default=NULL,nargs="+",
                   help="Subset of clusters to use for conversion [default: all clusters]")
parser$add_argument("-d", "--metadata_idents", type="character", default='seurat_clusters',
                   help="metadata column name to use for subseting")
parser$add_argument("-m", "--modalities", type="character", default=NULL,nargs = 2,
                   help="modalities to use for export - first is spliced, second unspliced")
parser$add_argument("-a", "--assay", type="character", default='GA',
                    help="assay to be used [requires common features within assay]")
parser$add_argument("-o", "--out", type="character", default='out.h5Seurat',
                    help="Path to the output .h5Seurat object [h5ad path determined from that]")
parser$add_argument("-t", "--pseudotime", type="character", default=NULL,
                    help="Path to seurat with pseudotime")


args <- parser$parse_args()

# Debug #########################
# args                 <- list()
# args$input           <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/multiple_modalities/ATAC_H3K27ac_H3K27me3/seurat_multimodal/peaks/Seurat_object.Rds'
# args$clusters        <- c('OPC','MOL')
# args$metadata_idents <- 'idents_short'
# args$modalities      <- c('ATAC','H3K27ac')
# args$assay           <- 'GA'
# args$out             <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/multiple_modalities/ATAC_H3K27ac_H3K27me3/seurat_multimodal/peaks/scvelo/ATAC_H3K27ac/Seurat_objet.h5seurat'
# args$pseudotime      <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/multimodal_data/WNN/seurat/pseudotime/Seurat_object_WNN_pseudotime.Rds'
###################################

seurat.ls <- readRDS(file=args$input)

# Remove merged object
if('merged' %in% names(seurat.ls)){
  seurat.ls <- seurat.ls[-match('merged',names(seurat.ls))]
}

# Filter common cells accross all modalities only
common.cells <- lapply(seurat.ls,Cells)

# If specified in arguments, filter specific clusters only 
if(!is.null(args$clusters)){
  common.cells <- lapply(seurat.ls,function(x){
    cells         <- x@meta.data[,args$metadata_idents]
    names(cells)  <- rownames(x@meta.data)
    cells         <- cells[cells %in% args$clusters]
    return(names(cells))
  })
}

# Intersect cells from various modalities
common.cells <- purrr::reduce(common.cells,intersect)

# Create new seurat object
bm              <-  CreateSeuratObject(seurat.ls[[args$modalities[1]]][[args$assay]][,common.cells],
                                       meta.data = seurat.ls[[args$modalities[1]]]@meta.data,assay='unspliced')
bm[['spliced']]  <- CreateAssayObject(seurat.ls[[args$modalities[2]]][[args$assay]]@counts[,common.cells],)

# Add pseudotime
if(!is.null(args$pseudotime)){
  seurat.pt <- readRDS(args$pseudotime)
  bm        <- AddMetaData(bm,metadata = seurat.pt$pt,col.name='pt')
}

# Run dimensionality reduction
bm <- SCTransform(object = bm, assay = "spliced")
bm <- RunPCA(object = bm, verbose = FALSE)

bm <- FindNeighbors(object = bm, dims = 1:20)
bm <- FindClusters(object = bm)
bm <- RunUMAP(object = bm, dims = 1:20)
DimPlot(bm,label=TRUE,group.by='idents_short') + NoLegend()

common.genes <- lapply(bm@assays, rownames)
common.genes <- purrr::reduce(common.genes,intersect)
common.genes <- intersect(common.genes, rownames(GetAssayData(bm,slot = 'scale.data',assay='SCT')))

bm <- bm[common.genes,]
bm@meta.data[,args$metadata_idents] <- as.character(bm@meta.data[,args$metadata_idents])

saveRDS(bm,file = paste0(args$out,'.Rds'))
SaveH5Seurat(bm, filename = args$out,overwrite=T)
Convert(args$out, dest = "h5ad",overwrite=T)

