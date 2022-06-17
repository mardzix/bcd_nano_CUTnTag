library(Seurat)
library(Signac)
library(argparse)
library(slingshot)
library(ggplot2)
library(viridis)

parser <- ArgumentParser()
parser$add_argument("--seurat", type="character",
                    help="Path to the seurat object")
parser$add_argument("--clusters", type="character",nargs="+",
                    help="Clusters to use for the pseudotime")
parser$add_argument("--reduction", type="character",default='umap',
                    help="Name of the dimensionality reduction")
parser$add_argument("--out", type="character",
                    help="Path to the seurat output")
args <- parser$parse_args()

# args <- list()
# args$seurat    <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/multimodal_data/WNN/seurat/Seurat_object_WNN.Rds'
# args$cluster   <- c('OPC','MOL')
# args$reduction <- 'wnn.umap'
# args$out <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/multimodal_data/WNN/seurat/pseudotime/Seurat_object_WNN_pseudotime.Rds'

seurat <- readRDS(file=args$seurat)
seurat.olg <- seurat[,seurat$idents_short %in% args$cluster]


# Remove outliers 
remove_outliers <- function(coords, nsd_cutoff = 5){
  nsd_vec <- (coords - mean(coords)) / sd(coords)
  nsd_vec <- nsd_vec[abs(nsd_vec) < nsd_cutoff]
  return(names(nsd_vec))
}

dim(seurat.olg)
DimPlot(seurat.olg,label = TRUE) + NoLegend()

seurat.olg <- seurat.olg[,remove_outliers(seurat.olg@reductions[[args$reduction]]@cell.embeddings[,1])]
seurat.olg <- seurat.olg[,remove_outliers(seurat.olg@reductions[[args$reduction]]@cell.embeddings[,1])]

dim(seurat.olg)
DimPlot(seurat.olg,label = TRUE) + NoLegend()

sshot <- slingshot(data = Embeddings(seurat.olg,reduction = args$reduction),
                   clusterLabels = seurat.olg$idents_short)

pt <- slingPseudotime(sshot)
seurat.olg$pt <- pt

p1 <- FeaturePlot(seurat.olg,features = 'pt') + scale_color_viridis_c()
ggsave(filename = paste0(args$out,'.png'),width=6,height=6)

saveRDS(object = seurat.olg,file = args$out)