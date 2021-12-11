library(Seurat)
library(Signac)
library(argparse)
library(rtracklayer)

# Tested for H3K27ac-H3K27me3 combination
# Would not work for ATAC probably

########### Arguments parser
cat("*** Parsing arguments\n")

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo',
                    help="seurat object")
parser$add_argument("-m", "--modalities", type="character", default=c('H3K27ac','H3K27me3'),
                    help="mod1",nargs="+")
parser$add_argument("-b", "--bed", type="character", default='foo',
                    help="path to the folder with bed files to be used")
parser$add_argument("-o", "--output_folder", type="character", default='foo',
                    help="output folder for files")
args <- parser$parse_args()

###########################
# 
args <- list()
args$input         <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/single_modality/H3K27ac/seurat_5000/Seurat_object_clustered_renamed.Rds'
args$modalities    <- c("H3K27ac","H3K27me3")
args$bed           <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/nbiotech_data/signal_matrix/top_peaks/'
args$output_folder <- '/data/proj/GCB_MB/bcd_CT/single-cell/results/'

fragments.gr <- lapply(args$modalities,function(x){
  fragments <- paste0('/data/proj/GCB_MB/bcd_CT/single-cell/results/single_modality/',x,'/fragments/fragments.tsv.gz')
  fragments
})
names(fragments.gr) <- args$modalities
##################
# Merged file
fragments.gr$merged <- paste0('/data/proj/GCB_MB/bcd_CT/single-cell/results/multiple_modalities/',paste0(args$modalities,collapse = '_'),'/fragments_merged_for_QC/fragments.tsv.gz')

################ -> Continue from here 
seurat    <- readRDS(file=args$input)

bed.ls       <- lapply(args$modalities,function(x){
  rtracklayer::import(paste0('/data/proj/GCB_MB/bcd_CT/single-cell/results/nbiotech_data/signal_matrix/top_peaks/peaks_',x,'_fragments.bed'))
  })


scores.ls <- lapply(bed.ls,function(bed_file){
  bed.matrix.ls <- lapply(fragments.gr,function(fragments_file){
    fragments.object <- CreateFragmentObject(path = fragments_file,cells = colnames(seurat))
    bed.matrix <- FeatureMatrix(fragments = fragments.object,
                                cells = colnames(seurat),
                                features = bed_file)
    bed.metadata <- colSums(bed.matrix)
    bed.metadata
  })
  bed.matrix.ls.norm <- lapply(bed.matrix.ls[1:length(bed.matrix.ls)-1],function(x){
    x/bed.matrix.ls[[length(bed.matrix.ls)]]
  })
  bed.matrix.ls.norm
})


scores.df <- lapply(scores.ls,function(x){
  n <- names(x)
  x <- purrr::reduce(.x = x,.f = cbind)
  colnames(x) <- n
  x
})
names(scores.df) <- args$modalities

for(x in names(scores.df)){
  seurat <- AddMetaData(seurat,metadata = scores.df[[x]],col.name = paste0(x,'_',colnames(scores.df[[x]])))
}


seurat <- AddMetaData(object = seurat,
                      metadata = bed.metadata[colnames(seurat)]/10^(seurat$logUMI),
                      col.name = basename(args$bed))

FeaturePlot(seurat,basename(args$bed)) + scale_color_viridis_c()



########################################################################################
n <- c('single_H3K27ac','single_H3K27me3','multiple_H3K27ac','multiple_H3K27me3')

fragments <- list('/data/proj/GCB_MB/bcd_CT/single-cell/results/nbiotech_data/fragments/H3K27ac_fragments.tsv.gz',
                  '/data/proj/GCB_MB/bcd_CT/single-cell/results/nbiotech_data/fragments/H3K27me3_fragments.tsv.gz',
                  '/data/proj/GCB_MB/bcd_CT/single-cell/results/single_modality/H3K27ac/fragments/fragments.tsv.gz',
                  '/data/proj/GCB_MB/bcd_CT/single-cell/results/single_modality/H3K27me3/fragments/fragments.tsv.gz')

seurat   <- list('/data/proj/GCB_MB/bcd_CT/single-cell/results/nbiotech_data/data/seurat/H3K27ac_seurat_object.Rds',
                 '/data/proj/GCB_MB/bcd_CT/single-cell/results/nbiotech_data/data/seurat/H3K27me3_seurat_object.Rds',
                 '/data/proj/GCB_MB/bcd_CT/single-cell/results/single_modality/H3K27ac/seurat_5000/Seurat_object_clustered_renamed.Rds',
                 '/data/proj/GCB_MB/bcd_CT/single-cell/results/single_modality/H3K27me3/seurat_5000/Seurat_object_clustered_renamed.Rds')

names(fragments) <- n
names(seurat)    <- n

fragments.ls <- lapply(fragments,function(x){
  rtracklayer::import(x,format='bed')
})

names(fragments.ls) <- n

bed.ls       <- lapply(args$modalities,function(x){ # TODO FIX
  rtracklayer::import(paste0('/data/proj/GCB_MB/bcd_CT/single-cell/results/nbiotech_data/signal_matrix/top_peaks/peaks_',x,'_all_fragments.bed'))
})
names(bed.ls) <- args$modalities


seurat.ls <- lapply(seurat,readRDS)
names(seurat.ls)<-n


cells.ls  <- lapply(seurat.ls,colnames)
names(cells.ls) <- n

cells.all <- unique(purrr::reduce(cells.ls,c))

fragments.ls.valid <- lapply(fragments.ls,function(x){
  x[x$name %in% cells.all]
})




findOverlapsGR <- function(bed,fragments,cells){
  op               <- findOverlaps(query = fragments,subject = bed)
  op.query         <- unique(queryHits(op))
  fragments.op     <- fragments[op.query]
  cells.op         <- table(fragments.op$name)
  cells.all        <- table(fragments$name)
  
  cells.not.in.op <- setdiff(union(names(cells.op),names(cells.all)),intersect(names(cells.op),names(cells.all)))
  cells.op[cells.not.in.op] <- 0
  return(cells.op/cells.all)
}

cells.stats <- list()

for(p in names(bed.ls)){
  cells.stats[[p]] <- list()
  for(f in names(fragments.ls.valid)){
    cells.stats[[p]][[f]] <- findOverlapsGR(bed = bed.ls[[p]],
                                             fragments = fragments.ls.valid[[f]],
                                             cells = cells.ls[[f]] 
                                             )
    seurat.ls[[f]] <- AddMetaData(object = seurat.ls[[f]], cells.stats[[p]][[f]],col.name = paste0(p,'_fraction'))
  
  }
}

df.to.plot.ls <- lapply(names(seurat.ls),function(x){
  data.frame('sample' = x,
             'cluster' = seurat.ls[[x]]@active.ident, 
             'H3K27ac_fraction' = seurat.ls[[x]]$H3K27ac_fraction,
             'H3K27me3_fraction' = seurat.ls[[x]]$H3K27me3_fraction)
})


df.to.plot <- purrr::reduce(df.to.plot.ls,rbind)

ggplot(data=df.to.plot,aes(x=sample,fill=cluster,y=H3K27ac_fraction)) + geom_boxplot() + coord_cartesian(ylim=c(0,1))


