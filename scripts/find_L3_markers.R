library(Seurat)
library(Signac)
library(argparse)
library(dplyr)
library(purrr)
set.seed(1234)

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo',
                    help="seurat object")
parser$add_argument("-o", "--output", type="character", default='foo',
                    help="output csv file for all markers")

args <- parser$parse_args()

#######
# args<-list()
# args$input  <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/multimodal_data/single_modality/H3K27ac/seurat/peaks/Seurat_object_clustered_renamed.Rds'
# args$output <-'/data/proj/GCB_MB/bcd_CT/single-cell/results/multimodal_data/single_modality/H3K27ac/seurat/peaks/markers/idents_L3_in_L1_niches/markers.csv'

seurat   <- readRDS(args$input)


markers     <- list()
markers.top <- list()

DefaultAssay(seurat) <- 'peaks'

for(identity in levels(seurat$idents_L1)){
  seurat.small              <- seurat[,seurat$idents_L1==identity]
  seurat.small              <- SetIdent(object = seurat.small,cells = names(seurat.small$idents_L3),value = seurat.small$idents_L3)
  markers.small             <- FindAllMarkers(seurat.small)
  markers.small$L1_identity <- identity
  markers[[identity]]       <- markers.small
  markers.top[[identity]]   <- markers.small %>% group_by(cluster) %>% filter(avg_log2FC > 0) %>% top_n(n=50,wt = -p_val_adj)
}

markers.out      <- purrr::reduce(markers,rbind)
markers.top.out  <- purrr::reduce(markers.top,rbind)

write.csv(file=args$output)
write.csv(file=gsub(pattern = ".csv",replacement = "_top.csv",x = args$output))




