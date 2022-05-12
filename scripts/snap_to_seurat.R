library(SnapATAC)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(rtracklayer)
library(argparse)
set.seed(1234)

parser <- ArgumentParser()
parser$add_argument("-s", "--snap", type="character", default='foo', 
                    help="path to the snap objects",nargs="+")
parser$add_argument("-m", "--metadata", type="character", default='foo', 
                    help="path to the metadata tsv",nargs="+")
parser$add_argument("-n", "--ncells", type="double", default=50000, 
                    help="number of cells to get/fraction of cells to get (if <1)")
parser$add_argument("-o", "--output", type="character", default='Seurat.Rds', 
                    help="output file")
args <- parser$parse_args()


# Load metadata
snap.metadata <- read.csv(file=args$metadata,row.names=NULL,sep='\t')
rownames(snap.metadata) <- paste0(snap.metadata$CellID)


# Load snap file
snap <- createSnap(file = args$snap,sample = gsub(".snap","",basename(args$snap)))

# Find common cells between metadata and snap
snap.cells <- paste0(snap@sample,'.',snap@barcode)
meta.cells <- rownames(snap.metadata)
common.cells <- intersect(snap.cells,meta.cells)

# Sample fraction of the barcodes
if(args$ncells <=1){
  ncells <- length(common.cells)*args$ncells
} else {
  ncells <- args$ncells
}

if(ncells > length(snap@barcode) | (ncells > dim(snap.metadata)[1])){
  stop('too many cells to sample, reduce number of cells')
}

barcodes.pick <- sample(common.cells,ncells)
barcodes.pick.index <- which(paste0(snap@sample,'.',snap@barcode) %in% barcodes.pick)

# Subset the snap file
snap <- snap[barcodes.pick.index,]

# Add matrices
snap <- SnapATAC::addGmatToSnap(snap)
snap <- SnapATAC::addPmatToSnap(snap)

# Create chromatin assay
counts.matrix           <- t(snap@pmat)
colnames(counts.matrix) <- paste0(snap@sample,'.',snap@barcode)
chrom.assay <- CreateChromatinAssay(counts=counts.matrix,ranges=snap@peak)
seurat      <- CreateSeuratObject(counts = chrom.assay,project = 'bingren',assay = 'peaks',meta.data = snap.metadata)

# Export files
saveRDS(file=args$output,object = seurat)
saveRDS(file=gsub(x = args$output,pattern = '.Rds',replacement = '_snap.Rds'),object = snap)

# Export fragments
# Default SnapATAC::extractReads gives a segmentation foult, so this is a workaround that 

frags.gr.selected <- list()

for(file in unique(snap@file)){
  sample <- unique(snap@sample[snap@file ==file])
  barcode.list = as.character(tryCatch(barcode.list <- h5read(file, '/BD/name'), error = function(e) {print(paste("Warning @extractReadsFromOneCell: 'BD/name' not found in ",file)); return(vector(mode="character", length=0))}))
  pos.list = as.numeric(tryCatch(pos.list <- h5read(file, "FM/barcodePos"), error = function(e) {print(paste("Warning @readMetaData: 'FM/barcodePos' not found in ",file)); return(c())}))
  len.list = as.numeric(tryCatch(len.list <- h5read(file, "FM/barcodeLen"), error = function(e) {print(paste("Warning @readMetaData: 'FM/barcodeLen' not found in ",file)); return(c())}))
  
  barcodes <- snap@barcode[which(snap@file == file)]
  chrom    <- h5read(file=file, name = 'FM/fragChrom')
  start    <- h5read(file=file, name='FM/fragStart')
  lens     <- h5read(file=file, name='FM/fragLen')
  frags.gr <- GRanges(chrom, 
                      IRanges(start, start + lens - 1))
  
  for(bcd in barcodes){
    id = paste0(sample,'.',bcd)
    pos = pos.list[match(bcd, barcode.list)]
    len = len.list[match(bcd, barcode.list)]
    idx.arr <- seq(pos, pos + len - 1)
    frags.gr.selected[[id]]          <- frags.gr[idx.arr]
    frags.gr.selected[[id]]$barcode  <- paste0(sample,'.',bcd)
    # frags.gr.selected[[bcd]] <- sortSeqlevels(frags.gr.selected[[bcd]])
    # frags.gr.selected[[bcd]] <- sort(frags.gr.selected[[bcd]])
    # frags.gr.selected[[bcd]]$file <- file
  }
}

gr <- as(frags.gr.selected,'GRangesList')
gr <- unlist(gr)
gr <- sortSeqlevels(gr)
gr <- sort(gr)
rtracklayer::export(object = gr, con = gsub(x = args$output, pattern = '.Rds',replacement = '_fragments.bed'),format = 'bed')




