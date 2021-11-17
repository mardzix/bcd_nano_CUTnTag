library(argparse)
library(Seurat)
library(Signac)
library(rtracklayer)
library(ggplot2);
library(funr)


source(paste0(dirname(funr::sys.script()),"/func.R"))
ndim = 40

###############################
parser <- ArgumentParser()

parser$add_argument("--sample", type="character",
                    help="sample name [as in config file]")

parser$add_argument("--antibody", type="character",
                    help="antibody name [as in config file]")

parser$add_argument("--metadata", type="character",
                    help="path to the pre-processed metadata [output from pick_cells.R]")

parser$add_argument("--fragments", type="character",
                    help="path to the fragments.tsv.gz with matching cell_barcodes")

parser$add_argument("--peaks", type="character",
                    help="path to the peaks file")

parser$add_argument("-o", "--out_prefix", type="character", default="10000",
                    help="folder for the output files")

parser$add_argument("--genome_version", type="character",
                    help="Version of the genome [e.g. mm10]")

parser$add_argument("-w", "--window", type="integer", default="10000", 
                    help="width of a window")

args <- parser$parse_args()

saveRDS(args,'arguments.Rds')

genome_version <- args$genome_version
############################ Filter the dataset

metadata <- read.csv(file = args$metadata,stringsAsFactors = FALSE)
metadata <- metadata[metadata$passedMB,]
rownames(metadata) <- metadata$barcode

######## Create  fragments object
fragments.path <- args$fragments
fragments      <- CreateFragmentObject(path = fragments.path,
                                       cells = metadata$barcode,
                                       verbose = TRUE,validate.fragments = TRUE)
####################
# Load annotations #
####################

# Peaks
peaks_file = args$peaks
if (!file.exists(peaks_file)) {stop(paste0("Peaks file does not exist:: ", peaks_file))}

# peaks <- rtracklayer::import(peaks_file,format='bed')
peaks <- read.table(file='results/bcdCT_MB21_02/ATAC_TATAGCCT/peaks/SEACR/peaks.relaxed.bed',stringsAsFactors=FALSE)
peaks <- GRanges(seq=peaks$V1,ranges=IRanges(start=peaks$V2,end=peaks$V3),score=peaks$V4)

# Genes
genebodyandpromoter.coords.flat <- load_ensembl_annot(genome_version)
genes.key             <- genebodyandpromoter.coords.flat$name
names(genes.key)      <- GRangesToString(genebodyandpromoter.coords.flat)

# Promoters
promoter.coords <- load_ensembl_promoters(genome_version)
promoters.key <- promoter.coords$tx_name
names(promoters.key) <- GRangesToString(promoter.coords)

###################
# Create matrices #
###################

cat("*** Creating peaks matrix \n")
counts.matrix.peaks <- FeatureMatrix(fragments = fragments,
                                     features = peaks,
                                     cells = metadata$barcode)


cat("*** Creating genomic bin matrix \n")
counts.matrix.bins <- GenomeBinMatrix(fragments = fragments,
                                      genome = setNames(getChromInfoFromUCSC(genome_version)[,2],getChromInfoFromUCSC(genome_version)[,1]),
                                      binsize = args$window,
                                      cells = metadata$barcode)

cat("*** Creating gene activities matrix \n")
gene.matrix     <- FeatureMatrix(fragments = fragments,
                                 features = genebodyandpromoter.coords.flat ,
                                 cells = metadata$barcode)

rownames(gene.matrix) <- genes.key[rownames(gene.matrix)]
gene.matrix           <- gene.matrix[rownames(gene.matrix) != "",]

cat("*** Creating promoter activities matrix \n")
promoter.matrix <- FeatureMatrix(fragments = fragments,
                                 features = promoter.coords,
                                 cells = metadata$barcode)

rownames(promoter.matrix) <- promoters.key[rownames(promoter.matrix)]
promoter.matrix <- promoter.matrix[rownames(promoter.matrix) != "",]

########################## Create Seurat object
min_features = 1
min_cells    = 1

seurat_object <- CreateSeuratObject(counts = counts.matrix.bins,
                     project = args$sample,
                     assay = paste0('bins_',args$window),
                     meta.data = metadata,
                     min.features = min_features,
                     min.cells = min_cells)


seurat_object[['peaks']] <- CreateAssayObject(counts = counts.matrix.peaks[,colnames(counts.matrix.peaks) %in% colnames(seurat_object)])
seurat_object[['GA']]    <- CreateAssayObject(counts = gene.matrix[,colnames(gene.matrix) %in% colnames(seurat_object)])
seurat_object[['PA']]    <- CreateAssayObject(counts = promoter.matrix[,colnames(promoter.matrix) %in% colnames(seurat_object)])

# Filter blacklist cells
seurat_object$blacklist_ratio <- seurat_object$blacklist_region_fragments / seurat_object$all_unique_MB
seurat_object                 <- seurat_object[,seurat_object$blacklist_region_fragments < 5]

#new.metadata <- unlist(config$samples[[args$sample]])
#for (x in seq(new.metadata)) {
#  seurat_object <- AddMetaData(object = seurat_object,metadata = new.metadata[x],col.name = names(new.metadata)[x])
#}


######### Save the object
cat("*** Save the object \n")
dir.create(args$out_prefix,recursive = TRUE)
saveRDS(object = seurat_object,file = paste0(args$out_prefix,'Seurat_object.Rds'))

########################## Try Clustering using Seurat
cat("*** Clustering and dimensionality reduction \n")

DefaultAssay(seurat_object) <- paste0('bins_',args$window)
# DefaultAssay(seurat_object) <- "peaks"

# The dimensionality reduction might fail, especially for low complexity datasets, so put the code into try to still save the objects
try({
seurat_object <- RunTFIDF(seurat_object)
seurat_object <- FindTopFeatures(seurat_object,min.cutoff = 'q0')
seurat_object <- RunSVD(
  seurat_object,
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
  )


seurat_object <- RunUMAP(seurat_object, dims = 2:ndim, reduction = 'lsi')

seurat_object <- FindNeighbors(
  object = seurat_object,
  reduction = 'lsi',
  dims = 2:ndim
)


seurat_object <- FindClusters(
  object = seurat_object,
  #algorithm = "leiden",
  resolution = 0.3,
  verbose = TRUE
)


p1 <- DimPlot(seurat_object,group.by = 'ident')
p2 <- FeaturePlot(seurat_object,'blacklist_region_fragments')
p3 <- FeaturePlot(seurat_object,'logUMI')



ggsave(plot= p1 + p2 + p3,
       filename = paste0(args$out_prefix,'Seurat_clustering.png'),
       width=15,height=5)
})


######### Save the final object
cat("*** Save the object \n")
saveRDS(object = seurat_object,file = paste0(args$out_prefix,'Seurat_object.Rds'))

