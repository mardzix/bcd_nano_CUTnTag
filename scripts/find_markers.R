library(Seurat)
library(Signac)
library(argparse)
library(funr)

source(paste0(dirname(funr::sys.script()),"/func.R"))
set.seed(1234)




########### Arguments parser
cat("*** Parsing arguments\n")

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo',
                    help="seurat object")
parser$add_argument("-o", "--output_file", type="character", default='foo',
                    help="output folder")
parser$add_argument("-d", "--idents", type="character", default='active.ident',
                    help="identities to use for markers")
parser$add_argument("-g", "--genome", type="character", default='mm10',
                    help="genome version to be used [e.g. mm10]")
args <- parser$parse_args()

###################

cat("*** Loading seurat object \n")
seurat <- readRDS(args$input)
# seurat <- readRDS(file='/data/proj/GCB_MB/bcd_CT/single-cell/results/single_modality/ATAC/seurat_5000/Seurat_object_clustered_renamed.Rds')


if(args$idents != 'active.ident'){
  seurat <- SetIdent(object = seurat,cells = names(seurat@meta.data[,args$idents]),value = seurat@meta.data[,args$idents])
  }




cat("*** Finding markers \n")
markers <- FindAllMarkers(seurat)

genes        <- load_ensembl_annot(args$genome)
closest_gene <- ClosestFeature(object = seurat, regions = StringToGRanges(markers$gene),genes)

markers               <- markers[markers$gene %in% closest_gene$query_region,]
closest_gene          <- setNames(closest_gene$name,closest_gene$query_region)
markers$closest_gene  <- closest_gene[markers$gene]

# SetNames(closest_gene$name,closest_gene$query_region)


markers          <- markers[markers$p_val_adj < 0.05,]
markers.positive <- markers[markers$avg_log2FC > 0,  ]

cat("*** Export markers as csv \n")
write.csv(x = markers,file = args$output_file)
write.csv(x = markers,file = gsub(pattern = '\\.csv','_positive.csv',args$output_file))

