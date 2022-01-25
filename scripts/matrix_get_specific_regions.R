library(GenomicRanges)
library(argparse)
library(rtracklayer)
library(dplyr)

########### Arguments parser
cat("*** Parsing arguments\n")

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo',
                    help="matrix tsv from multiBigWigSummary")
parser$add_argument("-o", "--output_folder", type="character", default='foo',
                    help="output folder for files")
args <- parser$parse_args()

###########################
# Input file == .tsv file
# $1 = chromosome
# $2 = start
# $3 = end
# $4 = score_1
# $5 = score_2
##########################
# args$input <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/nbiotech_data/signal_matrix/matrix.tab'

quantile_cutoff <- 0.8
npeaks=10000

m <- read.table(file=args$input,comment.char = "@",header = TRUE)
m$index <- 1:dim(m)[1]
colnames(m) <- gsub(pattern = '.bw',replacement = '',x = colnames(m))

m$log_FCH <- log2(m[,4]/m[,5])
m.gr      <- GRanges(seqnames = m[,1],ranges = IRanges(start = m[,2],end = m[,3]))


m.file1 <- m[m[,4] > quantile(m[,4], quantile_cutoff) ,]
m.file2 <- m[m[,5] > quantile(m[,5], quantile_cutoff) ,]

m.file1 <- head(m.file1[order(m.file1$log_FCH,decreasing = TRUE),],npeaks)
m.file2 <- head(m.file2[order(m.file2$log_FCH,decreasing = FALSE),],npeaks)

m.gr.file1 <- m.gr[as.numeric(m.file1$index),]
m.gr.file2 <- m.gr[as.numeric(m.file2$index),]

# Export
dir.create(args$output_folder)
rtracklayer::export(m.gr.file1,con = paste0(args$output_folder,'/peaks_',colnames(m)[4],'.bed'),format='bed')
rtracklayer::export(m.gr.file2,con = paste0(args$output_folder,'/peaks_',colnames(m)[5],'.bed'),format='bed')

