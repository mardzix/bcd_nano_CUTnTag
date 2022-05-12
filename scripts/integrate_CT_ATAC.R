library(argparse)
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(GenomicRanges)
set.seed(1234)

parser <- ArgumentParser()
parser$add_argument("--reference", type="character",
                    help="Path to reference seurat object")
parser$add_argument("--reference_assay", type="character",default='GA',
                    help="Reference assay to use for integration")
parser$add_argument("--reference_fragments", type="character",default=NULL,
                    help="Path to the reference fragments file [if not already in seurat object]")
parser$add_argument("--query", type="character",
                    help="path to query seurat objecy")
parser$add_argument("--query_assay", type="character",default='GA',
                    help="Query assay to use for integration")
parser$add_argument("--query_fragments", type="character",default=NULL,
                    help="Path to the query fragments file [if not already in seurat object]")
parser$add_argument("--reference_group", type="character",
                    help="Reference group (metadata) to use for plotting")
parser$add_argument("--query_group", type="character",
                    help="Query group (metadata) to use for plotting")
parser$add_argument("--downsample_features", type="integer",default=50000,
                    help="Downsample features")
parser$add_argument("--out", type="character",
                    help="path to the output")
args <- parser$parse_args()


############
# args<-list()
# args$reference <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/scATAC_bingren/seurat/Seurat_ATAC_clustered.Rds'
# args$query     <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/multimodal_data/single_modality/ATAC/seurat/peaks/Seurat_object_clustered_renamed.Rds'
# args$reference_assay <- 'peaks'
# args$query_assay <- 'peaks'
# args$reference_fragments <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/scATAC_bingren/seurat/Seurat_ATAC_fragments.bed.gz'
# args$out <- 'integration_test.Rds'
# args$downsample_features = 5000
# args$query_group <- 'idents_short'
# args$reference_group <- 'MajorType'


#######
cat("*** Loading seurat objects \n")
seurat.query          <- readRDS(file=args$query)
seurat.reference      <- readRDS(file=args$reference)


cat("*** Using query assay",args$query_assay," and reference assay", args$reference_assay, " for integration\n")
DefaultAssay(seurat.query) <- args$query_assay
DefaultAssay(seurat.reference) <- args$reference_assay

seurat.query$integration_ident     <- 'query'
seurat.reference$integration_ident <- 'reference'

if (is(seurat.query[[args$query_assay]],'ChromatinAssay') && is(seurat.reference[[args$reference_assay]],'ChromatinAssay')){
  all.features.gr <- GenomicRanges::intersect(seurat.query[[args$query_assay]]@ranges, seurat.reference[[args$reference_assay]]@ranges)
  
  if(length(all.features.gr) > args$downsample_features){
    cat("*** Found",length(all.features.gr),"for integration -> too many, randomly sample ",args$downsample_features," features\n")
    all.features.gr <- sample(all.features.gr,args$downsample_features)
  }
  
  cat("*** Creating matrices for new set  features\n")

  if(length(Fragments(seurat.reference)) == 0){
    cat("*** Add new fragments file to the query",args$reference_fragments ,"\n")
    Fragments(seurat.reference) <- CreateFragmentObject(path=args$reference_fragments)
  }
  if(length(Fragments(seurat.query)) == 0){
    cat("*** Add new fragments file to the query",args$query_fragments ,"\n")
    Fragments(seurat.query) <- CreateFragmentObject(path=args$query_fragments)
  }
  matrix.integration.query <- FeatureMatrix(fragments = Fragments(seurat.query),
                                            features = all.features.gr,
                                            cells = Cells(seurat.query))
  matrix.integration.reference <- FeatureMatrix(fragments = Fragments(seurat.reference),
                                                features = all.features.gr,
                                                cells = Cells(seurat.reference))
  
  seurat.query.integrated     <- CreateSeuratObject(counts = matrix.integration.query,assay = 'integration',meta.data = seurat.query@meta.data,min.cells=10,min.features=10)
  seurat.reference.integrated <- CreateSeuratObject(counts = matrix.integration.reference, assay = 'integration',meta.data = seurat.reference@meta.data,min.cells=10,min.features=10)
  
  DefaultAssay(seurat.reference.integrated) <- 'integration'
  DefaultAssay(seurat.query.integrated) <- 'integration'
  
  seurat.query.integrated     <- seurat.query.integrated %>% FindTopFeatures() %>% RunTFIDF()
  seurat.reference.integrated <- seurat.reference.integrated %>% FindTopFeatures() %>% RunTFIDF()
} 

all.features    <- intersect(rownames(seurat.reference.integrated),rownames(seurat.query.integrated))
all.features.gr <- StringToGRanges(all.features)


# Create transfer anchors using all the features
transfer.anchors <- FindTransferAnchors(
  reference = seurat.reference.integrated,
  query = seurat.query.integrated,
  reduction = 'cca',
  k.filter = 200,features = all.features
)

refdata <- GetAssayData(seurat.reference.integrated, slot = "data")
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = seurat.query[["lsi"]],dims = 2:50)

seurat.query.integrated[['integrated']]     <- imputation
seurat.reference.integrated[['integrated']] <- GetAssay(seurat.reference.integrated,DefaultAssay(seurat.reference.integrated))

if("GA" %in% Assays(seurat.query)){seurat.query.integrated[['GA']] <- seurat.query[['GA']]}
if("GA" %in% Assays(seurat.reference)){seurat.reference.integrated[['GA']] <- seurat.reference.integrated[['GA']]}

coembed               <- merge(x = seurat.reference.integrated, y = seurat.query.integrated)
DefaultAssay(coembed) <- 'integrated'

coembed <- RunSVD(coembed, features = all.features, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 2:40,reduction='lsi')

p1 <- DimPlot(coembed[,coembed$integration_ident=='reference'],label=TRUE,group.by = args$reference_group) + NoLegend()
p2 <- DimPlot(coembed[,coembed$integration_ident=='query'],label=TRUE,group.by=args$query_group) + NoLegend()
p1+p2

ggsave(filename = paste0(args$out,'.png'),plot = p1+p2,width=20,height=10)
saveRDS(object = coembed, file = args$out)
