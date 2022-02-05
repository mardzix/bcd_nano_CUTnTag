library(argparse)
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)

parser <- ArgumentParser()
parser$add_argument("--reference", type="character",
                    help="Path to reference seurat object")
parser$add_argument("--reference_assay", type="character",default='GA',
                    help="Reference assay to use for integration")
parser$add_argument("--query", type="character",
                    help="path to query seurat objecy")
parser$add_argument("--query_assay", type="character",default='GA',
                    help="Query assay to use for integration")
parser$add_argument("--out", type="character",
                    help="path to the output")
args <- parser$parse_args()

#######
cat("*** Loading seurat objects \n")
seurat.query          <- readRDS(file=args$query)
seurat.reference      <- readRDS(file=args$reference)

# args<-list()
# args$reference <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/scATAC_bingren/seurat/Seurat_merged_clustered.Rds'
# args$query     <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/multimodal_data/single_modality/ATAC/seurat/peaks/Seurat_object_clustered_renamed.Rds'
# args$reference_assay <- 'GA'
# args$query_assay <- 'GA'

cat("*** Using query assay",args$query_assay," and reference assay", args$reference_assay, " for integration\n")
DefaultAssay(seurat.query) <- args$query_assay
DefaultAssay(seurat.reference) <- args$reference_assay

seurat.query$integration_ident     <- 'query'
seurat.reference$integration_ident <- 'reference'


common.genes <- intersect(rownames(seurat.reference),rownames(seurat.query))
cat("*** Found",length(common.genes),"that will be used for integration\n")

transfer.anchors <- FindTransferAnchors(
  reference = seurat.reference,
  query = seurat.query,
  reduction = 'cca',
  k.filter = 200,features = common.genes
)


genes.use <- VariableFeatures(seurat.reference)
refdata <- GetAssayData(seurat.reference, assay = args$reference_assay, slot = "data")[genes.use, ]

imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = seurat.query[["lsi"]],dims = 2:40)

seurat.query[['integrated']]     <- imputation
seurat.reference[['integrated']] <- seurat.reference[[args$reference_assay]]


coembed               <- merge(x = seurat.reference, y = seurat.query)
DefaultAssay(coembed) <- 'integrated'

coembed <- RunSVD(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 2:40,reduction='lsi')

p1 <- DimPlot(coembed[,coembed$integration_ident=='reference'],label=TRUE) + NoLegend()
p2 <- DimPlot(coembed[,coembed$integration_ident=='query'],label=TRUE) + NoLegend()
p1+p2

ggsave(filename = paste0(args$out,'.png'),plot = p1+p2,width=20,height=10)
saveRDS(object = coembed, file=args$out)
