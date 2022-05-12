library(argparse)
library(Seurat)
library(Signac)
library(rtracklayer)
library(ggplot2);
library(funr)


source(paste0(dirname(funr::sys.script()),"/func.R"))
# ndim = 40

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

parser$add_argument("--peaks", type="character",default=FALSE,
                    help="path to the peaks file")

parser$add_argument("-w", "--window", type="integer", default="10000",
                    help="width of a window for the BINSxCELLS matrix (if no peaks file is supplied)")

parser$add_argument("-o", "--out_prefix", type="character", default="10000",
                    help="folder for the output files")

parser$add_argument("--genome_version", type="character",
                    help="Version of the genome for gene activity matrix [e.g. mm10]")

parser$add_argument("-n", "--ndim", type="integer", default="30", 
                    help="number of LSI dimensions to use for initial clustering")


args <- parser$parse_args()

saveRDS(args,'arguments.Rds')

######################## Assay
if(is.character(args$peaks)){
    assay = 'peaks'
} else {
    assay = paste0('bin_',args$window)
}

############################ Filter the dataset

metadata <- read.csv(file = args$metadata,stringsAsFactors = FALSE)
metadata <- metadata[metadata$passedMB,]
rownames(metadata) <- metadata$barcode

######## Create  fragments object
fragments.path <- args$fragments
fragments      <- CreateFragmentObject(path = fragments.path,
                                       cells = metadata$barcode,
                                       verbose = TRUE,validate.fragments = TRUE)
#########################
# Load gene annotations #
#########################

# Peaks
if(assay == 'peaks'){
    if (!file.exists(args$peaks)) {stop(paste0("Peaks file does not exist:: ", peaks_file))}

    peaks <- read.table(file=args$peaks,stringsAsFactors=FALSE)
    peaks <- GRanges(seq=peaks$V1,ranges=IRanges(start=peaks$V2,end=peaks$V3),score=peaks$V4)


    cat("*** Creating peaks matrix \n")
    counts.matrix <- FeatureMatrix(fragments = fragments,
                                   features = peaks,
                                   cells = metadata$barcode)

    }

# Genomic bin matrix
if(assay != 'peaks'){
    cat("*** Creating genomic bin matrix \n")
    counts.matrix <- GenomeBinMatrix(fragments = fragments,
                                     genome = setNames(getChromInfoFromUCSC(args$genome_version)[,2],getChromInfoFromUCSC(args$genome_version)[,1]),
                                     binsize = args$window,
                                     cells = metadata$barcode)
    }



# Gene activity matrix
genebodyandpromoter.coords.flat        <- load_ensembl_annot(args$genome_version)
names(genebodyandpromoter.coords.flat) <- genebodyandpromoter.coords.flat$name

genes.key             <- genebodyandpromoter.coords.flat$name
names(genes.key)      <- GRangesToString(genebodyandpromoter.coords.flat)

cat(paste0("*** Creating gene activities matrix for genome ",args$genome_version," \n"))
gene.matrix     <- FeatureMatrix(fragments = fragments,
                                 features = genebodyandpromoter.coords.flat ,
                                 cells = metadata$barcode)

rownames(gene.matrix) <- genes.key[rownames(gene.matrix)]
gene.matrix           <- gene.matrix[rownames(gene.matrix) != "",]

gene.matrix          <- gene.matrix[,colSums(gene.matrix) > 0]


########################## Create Seurat object
min_features = 1
min_cells    = 1

chromatin.assay <- CreateChromatinAssay(counts = counts.matrix[,colnames(gene.matrix)],
                                        min.cells = min_cells,
                                        min.features = min_features,
                                        fragments = fragments.path,
                                        genome = args$genome_version)

seurat_object   <- CreateSeuratObject(counts= chromatin.assay,
                                      project = args$sample,
                                      assay = assay,
                                      meta.data = metadata)


###### Add some metadata
seurat_object$modality <- args$antibody
seurat_object$sample   <- args$sample

############### Add GA assay to the object
seurat_object[['GA']] <- CreateAssayObject(counts = gene.matrix[,Cells(seurat_object)],
                                           min.cells = min_cells,
                                           min.features = min_features)

################### Remove cells with very high number of reads
cat("*** Removing cells with very high UMI count \n")
logUMI_cutoff_high <- quantile(seurat_object$logUMI, 0.95)
logUMI_cutoff_low  <- quantile(seurat_object$logUMI, 0.05)

cat("*** Removing ",sum(seurat_object$logUMI < logUMI_cutoff_low)," cells for having UMI less than ",logUMI_cutoff_low,", \n")
cat("*** Removing ",sum(seurat_object$logUMI > logUMI_cutoff_high)," cells for having UMI more than ",logUMI_cutoff_high,", \n")

seurat_object <- seurat_object[,seurat_object$logUMI > logUMI_cutoff_low]
seurat_object <- seurat_object[,seurat_object$logUMI < logUMI_cutoff_high]

######### Save the object
cat("*** Save the object \n")
dir.create(args$out_prefix,recursive = TRUE)
saveRDS(object = seurat_object,file = paste0(args$out_prefix,'Seurat_object.Rds'))

########################## Try Clustering using Seurat
cat("*** Clustering and dimensionality reduction \n")


# The dimensionality reduction might fail, especially for low complexity datasets, so put the code into try to still save the objects

# Get assays and reorder
assays_all <- names(seurat_object@assays)
assays_all <- assays_all[order(c(grepl('bin|peak',assays_all)),decreasing = FALSE)]


for(assay in assays_all) {
  try({
    cat("*** Clustering on assay ",assay," \n")
    DefaultAssay(seurat_object) <- assay
    seurat_object <- RunTFIDF(seurat_object)
    seurat_object <- FindTopFeatures(seurat_object,min.cutoff = 'q0')
    seurat_object <- RunSVD(
      seurat_object,
      reduction.key = 'LSI_',
      reduction.name = 'lsi',
      assay = assay
      )
    
    
    seurat_object <- RunUMAP(seurat_object, dims = 2:args$ndim, reduction = 'lsi',assay = assay)
    
    seurat_object <- FindNeighbors(
      object = seurat_object,
      reduction = 'lsi',
      dims = 2:args$ndim,
      assay = assay
    )
    
    
    seurat_object <- FindClusters(
      object = seurat_object,
      verbose = TRUE,
      assay=assay
    )
    
    
    p1 <- DimPlot(seurat_object,group.by = 'ident',label=TRUE) + NoLegend()
    p2 <- FeaturePlot(seurat_object,'logUMI')
    
    
    
    ggsave(plot= p1 + p2,
           filename = paste0(args$out_prefix,'Seurat_clustering_',assay,'.png'),
           width=10,height=5)
    })
}

######### Save the final object
cat("*** Save the object \n")
saveRDS(object = seurat_object,file = paste0(args$out_prefix,'Seurat_object.Rds'))

