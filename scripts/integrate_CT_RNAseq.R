cat("*** Loading libraries\n")
library(Seurat)
library(Signac)
library(ggplot2)
library(argparse)


########### Arguments parser
cat("*** Parsing arguments\n")

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo',
                    help="multimodal seurat object to integrate")

parser$add_argument("-r", "--rna_seq", type="character", default='foo',
                    help="scrna_seq dataset to use")

parser$add_argument("-a", "--assay", type="character", default='GA',
                    help="Assay to use for integration [CT data]")

parser$add_argument("-m", "--modality", type="character", default=c("H3K27ac","ATAC","H3K4me3","H3K27me3"),
                    help="Modality to use for integration [CT data]")

parser$add_argument("-o", "--output", type="character", default='foo',
                    help="output file")
args <- parser$parse_args()

saveRDS(object = args,file = 'arguments.Rds')

cat("*** Loading data \n")
# rna.seq <- readRDS(file='/data/proj/GCB_MB/single-cell-CUT-Tag/nbiotech_paper/analysis/results/Sten_RNA/clustering/01.clustering_20000cells.Rds')
# seurat  <- readRDS(file='/data/proj/GCB_MB/bcd_CT/single-cell/results/multiple_modalities/H3K27ac_H3K27me3/seurat_5000/Seurat_object_clustered.Rds')

rna.seq <- readRDS(file=args$rna_seq)
seurat  <- readRDS(file=args$input)

# Convert single modality object to a named list
if(length(seurat) == 1){
  seurat <- setNames(object = list(seurat),nm = unique(seurat$modality))
}

# Add experiment metadata
cat("*** Add experiment metadata \n")
rna.seq$experiment <- "scRNAseq"
seurat$experiment  <- "bcdCT"

# Remove PNS cells
cat("*** Remove PNS cells \n")
rna.seq <- rna.seq[,-grep("Enteric",rna.seq$TaxonomyRank4)]
rna.seq <- rna.seq[,-grep("Peripheral",rna.seq$TaxonomyRank4)]

# Find assay for integration
integration_modality <- args$modality[min(which(args$modality %in% names(seurat)))]
cat("*** Modality to use for integration:",integration_modality,"\n")

integration_assay    <- args$assay
cat("*** Assay to use for integration:",integration_assay,"\n")

seurat.to.integrate               <- seurat[[integration_modality]]
DefaultAssay(seurat.to.integrate) <- integration_assay
seurat.to.integrate$experiment  <- "bcdCT"

cat("*** Features in CT data\n")
head(rownames(seurat.to.integrate))
cat("*** Features in scrna_seq data \n")
head(rownames(rna.seq))

common.features <- table(c(rownames(seurat.to.integrate),rownames(rna.seq)))
common.features <- names(common.features[common.features==2])

cat("*** Common features\n")
head(common.features)
length(common.features)


# Integrate
cat("*** Starting integration\n")
transfer.anchors <- FindTransferAnchors(
  reference = rna.seq,
  query = seurat.to.integrate,
  reduction = 'cca',
  query.assay = integration_assay,reference.assay = 'RNA',
  features = common.features
)


genes.use <- VariableFeatures(rna.seq)
refdata <- GetAssayData(rna.seq, assay = "RNA", slot = "data")[genes.use, ]

cat("*** Imputing data\n")
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = seurat.to.integrate[["lsi"]],dims = 2:50)

seurat.to.integrate[['RNA']] <- imputation
cat("*** Merging objects\n")
coembed <- merge(x =rna.seq , y = seurat.to.integrate)

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)


p1 <- DimPlot(seurat.to.integrate,pt.size=0.4,repel=TRUE,label=TRUE) + NoLegend() + ggtitle('multimodal CUT&Tag un-integrated')
p2 <- DimPlot(rna.seq,label=TRUE,repel=TRUE,pt.size=0.2,group.by="TaxonomyRank4") + NoLegend() + ggtitle('mouse brain atlas scrna_seq un-integrated')

p3 <- DimPlot(coembed[,coembed$orig.ident != 'bcdCT'],group.by='TaxonomyRank3',label=TRUE,pt.size = 0.2,repel = TRUE,) + NoLegend()
p4 <- DimPlot(coembed[,coembed$orig.ident == 'bcdCT'],label=TRUE,pt.size=0.2) + NoLegend()
p5 <- DimPlot(coembed[,coembed$orig.ident != 'bcdCT'],group.by='TaxonomyRank4',label=TRUE,pt.size = 0.2,repel=TRUE) + NoLegend()
p6 <- DimPlot(coembed[,coembed$orig.ident != 'bcdCT'],group.by='ClusterName',label=TRUE,pt.size = 0.2,repel=TRUE) + NoLegend()

pdf(file = paste0(dirname(args$output),'/integration_scRNAseq.pdf'),width=20,height=10)
p1+p2
p3 + p5 + p4
p4+p6
DimPlot(coembed,group.by='experiment')
dev.off()

dir.create(dirname(args$output),recursive = TRUE)
saveRDS(object = coembed,file = args$output)