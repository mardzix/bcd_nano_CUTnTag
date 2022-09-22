library(argparse)
library(GenomicRanges)
library(rtracklayer)

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo',
                    help="matrix output from multiBigwigSummary [deeptools]")
parser$add_argument("-o", "--output", type="character", default='foo',
                    help="path to output bed file")
parser$add_argument("-c", "--column", type="integer", default=1,
                    help="which signal column to use [as defined from file order in compute matrix]")
parser$add_argument("-q", "--quantile", type="double", default='0.2',
                    help="quantile fraction to use for cutoff")
parser$add_argument("-l", "--less", type="logical", default=FALSE,
                    help="less than quantile")
parser$add_argument("-g", "--greater", type="logical", default=FALSE,
                    help="greater than quantile")
args <- parser$parse_args()


m     <- read.table(file=args$input)

m.idx           <- m[,args$column + 3]
quantile.cutoff <- quantile(m.idx,args$quantile)

if(args$less){m.idx <- which(m.idx < quantile.cutoff)}
if(args$greater){m.idx <- which(m.idx > quantile.cutoff)}

m.bed <- m[m.idx,]
m.bed <- GRanges(seqnames = m.bed[,1], ranges = IRanges(start = m.bed[,2], end = m.bed[,3]))

rtracklayer::export(object = m.bed, con = args$output)