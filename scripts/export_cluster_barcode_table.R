library(Seurat)
library(Signac)
library(argparse)
set.seed(1234)

source(paste0(dirname(funr::sys.script()),"/func.R"))

########### Arguments parser
cat("*** Parsing arguments\n")

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo',
                    help="seurat object")
parser$add_argument("-o", "--output", type="character", default='foo',
                    help="output csv file")
parser$add_argument("-d", "--idents", type="character", default='active.ident',
                    help="identities to use for table")

args <- parser$parse_args()

###################
seurat <- readRDS(file=args$input)

if(args$idents != 'active.ident'){
  seurat <- SetIdent(object = seurat,cells = names(seurat@meta.data[,args$idents]),value = seurat@meta.data[,args$idents])
  }


cluster.annotations           <- data.frame(cluster=seurat@active.ident,
                                            barcode=names(seurat@active.ident))
colnames(cluster.annotations) <- paste0("#",colnames(cluster.annotations))
write.csv(x = cluster.annotations, file=args$output,row.names = FALSE,quote = FALSE)
