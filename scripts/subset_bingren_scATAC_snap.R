library(SnapATAC)
library(purrr)
set.seed(1234)


########### Arguments parser
cat("*** Parsing arguments\n")

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo',
                    help="seurat object")
parser$add_argument("-o", "--output", type="character", default='foo',
                    help="output csv file")
parser$add_argument("-f", "--fraction", type="numeric", default=0.05,
                    help="Fraction of cells  to extract")
args <- parser$parse_args()

cell.ids.ls <- lapply(args$input,function(x){
  cat("*** Reading barcodes from file",basename(x),"\n")
  snap <- SnapATAC::createSnap(file=x, sample=gsub("\\.snap","",basename(x)))
  return(snap@barcode)
})

all.cells  <- unlist(cell.ids.ls)
cells.pick <- sample(all.cells,length(all.cells)*args$fraction)

cell.ids.ls <- lapply(args$input,function(x){
  cat("*** Loading cells from snap file",basename(x),"\n")
  snap <- SnapATAC::createSnap(file=x, sample=gsub("\\.snap","",basename(x)))
  return(snap@barcode)
})


