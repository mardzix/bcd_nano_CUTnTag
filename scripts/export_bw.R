library(Seurat)
library(Signac)
library(argparse)
library(rtracklayer)
library(GenomeInfoDb)

set.seed(1234)

source(paste0(dirname(funr::sys.script()),"/func.R"))

########### Arguments parser
cat("*** Parsing arguments\n")

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo',
                    help="seurat object")
parser$add_argument("-f", "--fragments", type="character", default='foo',
                    help="path to the corresponding fragments file")
parser$add_argument("-o", "--output_folder", type="character", default='foo',
                    help="output folder")
parser$add_argument("-d", "--idents", type="character", default='active.ident',
                    help="identities to use for export")

args <- parser$parse_args()

###################
# args <- list()
# args$input <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/single_modality/H3K27ac/seurat_5000/Seurat_object_clustered_renamed.Rds'
# saveRDS(args,'arguments.Rds')
###################

cat("*** Loading data \n")
seurat_object <- readRDS(file=args$input)

chrom.sizes <- getChromInfoFromUCSC('mm10',assembled.molecules.only=TRUE)
fragments <- read.delim(gzfile(args$fragments),sep='\t',row.names=NULL,header=FALSE)
fragments <- GRanges(seqnames=fragments$V1, ranges = IRanges(start=fragments$V2,end=fragments$V3),name=fragments$V4, score=fragments$V5)

dir.create(args$output_folder,recursive = TRUE)

if(args$idents != 'active.ident'){
  seurat_object@active.ident <- seurat_object@meta.data[,args$idents]
}

clusters_to_use <- levels(seurat_object@active.ident)

lapply(clusters_to_use,function(x){
  exportBW(object = seurat_object,
           cluster = x,
           fragments = fragments,
           path = paste0(args$output_folder,'/cluster_',x,'.bw'),chrom.sizes=chrom.sizes)
})