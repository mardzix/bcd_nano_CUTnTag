library(argparse)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(Signac)
set.seed(1234)

parser <- ArgumentParser()
parser$add_argument("--input", type="character",
                    help="Path to the markers csv file")
parser$add_argument("--nmarkers", type="integer",
                    help="Number of top marker regions to export")
parser$add_argument("--output", type="character",
                    help="Output bed file")
args <- parser$parse_args()

markers <- read.csv(file=args$input)
markers <- markers %>% group_by(cluster) %>% top_n(n = args$nmarkers,wt = -log(p_val_adj))
markers.gr <- StringToGRanges(markers$gene)
markers.gr <- GenomicRanges::reduce(markers.gr)

rtracklayer::export(object = markers.gr,con = args$output,format = 'bed')