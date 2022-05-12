library(argparse)
library(reticulate)
library(rtracklayer)
library(GenomicRanges)

parser <- ArgumentParser()
parser$add_argument("--matrix_txt", type="character",
                    help="Path to multiBigwigSummary txt matrix")
parser$add_argument("--matrix_npz", type="character",
                    help="Path to multiBigwigSummary npz matrix")
parser$add_argument("--output", type="character",
                    help="Path to outputs folder")
args <- parser$parse_args()



# args <- list()
# args$matrix_txt <- '/data/proj/GCB_MB/bcd_CT/single-cell/results//benchmarks/peaks/specificity_benchmark/peaks_matrix.txt'
# args$matrix_npz <- '/data/proj/GCB_MB/bcd_CT/single-cell/results//benchmarks/peaks/specificity_benchmark/peaks_matrix.npz'
# args$output <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/benchmarks/peaks/specificity_benchmark/peaks_specific/'

np <- reticulate::import('numpy')

m     <- read.table(file=args$matrix_txt)
m.np <- np$load(args$matrix_npz)

peaks <- GRanges(seqnames = m$V1,ranges = IRanges(start = m$V2,end = m$V3))
m <- m[,-1:-3]
colnames(m) <- m.np[['labels']]


m$log2_ratio <- log2(m$H3K27ac_Astrocytes.bw + 1e-4)/(m$H3K27me3_Astrocytes.bw + 1e-4)

K27ac_peaks  <- which(m$H3K27ac_Astrocytes.bw > quantile(m$H3K27ac_Astrocytes.bw,0.9) & m$log2_ratio > 4)
K27me3_peaks <- which(m$H3K27me3_Astrocytes.bw > quantile(m$H3K27me3_Astrocytes.bw,0.9) & m$log2_ratio < 4)

dir.create(args$output,recursive = TRUE)
rtracklayer::export(object = sort(peaks[K27ac_peaks]),con = paste0(args$output,'/peaks_K27ac.bed'))
rtracklayer::export(object = sort(peaks[K27me3_peaks]),con = paste0(args$output,'/peaks_K27me3.bed'))
