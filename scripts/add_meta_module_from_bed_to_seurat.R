library(Seurat)
library(Signac)
library(argparse)
library(rtracklayer)


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

  return(cells.op/cells.all[names(cells.op)])
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

ggplot(data=df.to.plot,aes(x=sample,y=H3K27ac_fraction)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim=c(0,0.3))


