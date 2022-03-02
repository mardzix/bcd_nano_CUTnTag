library(Seurat)
library(Signac)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(purrr)
library(slingshot)
library(cicero)
set.seed(1234)


seurat <- readRDS(file='/data/proj/GCB_MB/bcd_CT/single-cell/results/multimodal_data/WNN/seurat/Seurat_object_WNN.Rds')


seurat.olg <- seurat[,seurat$idents_short %in% c('OPC','MOL')]
seurat.olg <- seurat.olg[,Embeddings(seurat.olg, reduction = 'wnn.umap')[,1] < -7]

DimPlot(seurat,group.by='idents_short', label=TRUE) + NoLegend()
DimPlot(seurat.olg,group.by='idents_short', label=TRUE) + NoLegend()

sshot <- slingshot(data = Embeddings(seurat.olg,reduction = 'wnn.umap'),clusterLabels = seurat.olg$idents_short)
plot(sshot)
pt <- slingPseudotime(sshot)
seurat.olg$pt <- pt

FeaturePlot(seurat.olg,features = 'pt') + scale_color_viridis_c()

# imputed <- magic(data = seurat.olg)


modalities <- c("ATAC","H3K27ac","H3K27me3")

for(i in modalities){
  for(x in modalities){
    markers <- read.csv(file=paste0('/data/proj/GCB_MB/bcd_CT/single-cell/results/multimodal_data/single_modality/',i,'/seurat/peaks/markers/idents_short/markers_positive.csv'))
    markers.mol <- head(markers[markers$cluster == 'MOL',],100)
    
    DefaultAssay(seurat.olg) <- paste0('peaks_',x)
    op <- findOverlaps(query = StringToGRanges(rownames(seurat.olg)), StringToGRanges(markers.mol$gene))
    
    d <- GetAssayData(object = seurat.olg,slot = 'data',assay = paste0('peaks_',x))  
    MOL.score <- colSums(d[unique(queryHits(op)),])
    seurat.olg <- AddMetaData(object = seurat.olg,metadata = MOL.score,col.name = paste0('MOL.',x,'.score.',i))
  }
}


d.small <- d[unique(queryHits(op)),order(seurat.olg$pt)]


pheatmap(d.small,cluster_rows = FALSE,cluster_cols = FALSE,labels_col = FALSE, labels_row = NULL)


p1 <- ggplot(data = seurat.olg@meta.data,aes(x=pt,y=MOL.ATAC.score)) + geom_point() + geom_smooth(method='lm',formula = y ~ poly(x, 4))
p2 <- ggplot(data = seurat.olg@meta.data,aes(x=pt,y=MOL.H3K27ac.score)) + geom_point()+ geom_smooth(method='lm',formula = y ~ poly(x, 4))
p3 <- ggplot(data = seurat.olg@meta.data,aes(x=pt,y=MOL.H3K27me3.score)) + geom_point()+ geom_smooth(method='lm',formula = y ~ poly(x,4))


p1+p2+p3
############################################

heatmap.ls <- list()
modalities <- c("ATAC","H3K27ac","H3K27me3")

for(i in modalities){
  heatmap.ls[[i]] <- list()
  for(x in modalities){
    markers <- read.csv(file=paste0('/data/proj/GCB_MB/bcd_CT/single-cell/results/multimodal_data/single_modality/',i,'/seurat/peaks/markers/idents_short/markers.csv'))
    
    # markers.mol <- markers[markers$cluster %in% c('MOL','OPC'),]
    markers.mol <- markers[markers$cluster == 'MOL',]
    markers.mol <- markers.mol[markers.mol$avg_log2FC > 0,]
    markers.mol <- markers.mol[order(markers.mol$avg_log2FC,decreasing=TRUE),]
    markers.mol <- markers.mol[markers.mol$p_val_adj < 0.05 ,]
    markers.mol <- head(markers.mol,200)
    
    DefaultAssay(seurat.olg) <- paste0('peaks_',x)
    op <- findOverlaps(query = StringToGRanges(rownames(seurat.olg)), StringToGRanges(markers.mol$gene))
    
    d <- GetAssayData(object = seurat.olg,slot = 'data',assay = paste0('peaks_',x))  
    heatmap.ls[[i]][[x]] <- d[unique(queryHits(op)),order(seurat.olg$pt)]
  }
}

col_annot <- data.frame(row.names = rownames(seurat.olg@meta.data),
                        'cluster' = seurat.olg@meta.data[,'idents_short'],
                        'pt' = seurat.olg$pt)




pheatmap(heatmap.ls[['ATAC']][['ATAC']],cluster_cols = F,show_colnames = F, show_rownames =F, annotation_col = col_annot,treeheight_row = 0)
pheatmap(heatmap.ls[['ATAC']][['H3K27ac']],cluster_cols = F,show_colnames = F,show_rownames=F,annotation_col = col_annot,treeheight_row = 0)
pheatmap(heatmap.ls[['ATAC']][['H3K27me3']],cluster_cols = F,show_colnames = F,show_rownames=F,annotation_col = col_annot,treeheight_row = 0)



pheatmap(heatmap.ls[['H3K27ac']][['ATAC']],cluster_cols = F,show_colnames = F,show_rownames=F,annotation_col = col_annot)
pheatmap(heatmap.ls[['H3K27ac']][['H3K27ac']],cluster_cols = F,show_colnames = F,show_rownames=F,annotation_col = col_annot)
pheatmap(heatmap.ls[['H3K27ac']][['H3K27me3']],cluster_cols = F,show_colnames = F,show_rownames=F,annotation_col = col_annot)


pheatmap(heatmap.ls[['H3K27me3']][['ATAC']],cluster_cols = F,show_colnames = F,show_rownames=F,annotation_col = col_annot)
pheatmap(heatmap.ls[['H3K27me3']][['H3K27ac']],cluster_cols = F,show_colnames = F,show_rownames=F,annotation_col = col_annot)
pheatmap(heatmap.ls[['H3K27me3']][['H3K27me3']],cluster_cols = F,show_colnames = F,
         show_rownames=F,annotation_col = col_annot,scale = 'row')


#########################
a <- sort(kmeans(x=heatmap.ls[['H3K27ac']][['H3K27ac']],centers = 2)$cluster)
pheatmap(heatmap.ls[['H3K27ac']][['H3K27ac']][names(a),],cluster_cols = F,cluster_rows = F,show_colnames = F,show_rownames=F,annotation_col = col_annot)

df1 <-data.frame('score'           = colSums(heatmap.ls[['H3K27ac']][['H3K27ac']]), 
                 'pt'              = seurat.olg$pt[colnames(heatmap.ls[['H3K27ac']][['H3K27ac']])],
                 'modality'        = 'H3K27ac')

for(cluster in unique(a)){
  df1[,paste0('wave_',cluster)] <- colSums(heatmap.ls[['H3K27ac']][['H3K27ac']][names(a[a==cluster]),])
}

p1 <- ggplot(data=df1,aes(x=pt,y=score)) + geom_point() + geom_smooth(method='loess')
p1

library(reshape2)
df1.long <- reshape2::melt(df1,id = list('pt','modality','score'))
p1 <- ggplot(data=df1.long,aes(x=pt,y=value,col=variable)) + geom_point() + geom_smooth(method='loess')
p1



#######################################################################################
a <- sort(kmeans(x=heatmap.ls[['H3K27me3']][['H3K27me3']],centers = 2)$cluster)
pheatmap(heatmap.ls[['H3K27me3']][['H3K27me3']][names(a),],cluster_cols = F,cluster_rows = F,show_colnames = F,show_rownames=F,annotation_col = col_annot)
  #         scale='row',color = c(rep('white',10),colorRampPalette(c('white','red'))(40),rep('red',10)))

df1 <-data.frame('score'           = colSums(heatmap.ls[['H3K27me3']][['H3K27me3']]), 
                 'pt'              = seurat.olg$pt[colnames(heatmap.ls[['H3K27me3']][['H3K27me3']])],
                 'modality'        = 'H3K27ac')

for(cluster in unique(a)){
  df1[,paste0('wave_',cluster)] <- colSums(heatmap.ls[['H3K27me3']][['H3K27me3']][names(a[a==cluster]),])
}

p1 <- ggplot(data=df1,aes(x=pt,y=score)) + geom_point() + geom_smooth(method='loess')
p1

library(reshape2)
df1.long <- reshape2::melt(df1,id = list('pt','modality','score'))
p1 <- ggplot(data=df1.long,aes(x=pt,y=value,col=variable)) + geom_point() + geom_smooth(method='loess')
p1

#####################################
library(EnsDb.Mmusculus.v79)
library(ensembldb)


for(i in modalities){
  heatmap.ls[[i]] <- list()
  for(x in modalities){
    markers <- read.csv(file=paste0('/data/proj/GCB_MB/bcd_CT/single-cell/results/multimodal_data/single_modality/',i,'/seurat/peaks/markers/idents_short/markers.csv'))
    
    # markers.mol <- markers[markers$cluster %in% c('MOL','OPC'),]
    markers.mol <- markers[markers$cluster == 'MOL',]
    markers.mol <- markers.mol[markers.mol$avg_log2FC > 0,]
    markers.mol <- markers.mol[order(markers.mol$avg_log2FC,decreasing=TRUE),]
    markers.mol <- markers.mol[markers.mol$p_val_adj < 0.01 ,]
    markers.mol <- head(markers.mol,500)
    
    DefaultAssay(seurat.olg) <- paste0('peaks_',x)
    op <- findOverlaps(query = StringToGRanges(rownames(seurat.olg)), StringToGRanges(markers.mol$gene))
    
    d <- GetAssayData(object = seurat.olg,slot = 'data',assay = paste0('peaks_',x))  
    heatmap.ls[[i]][[x]] <- d[unique(queryHits(op)),order(seurat.olg$pt)]
  }
}

a <- sort(kmeans(x=heatmap.ls[['H3K27me3']][['H3K27me3']],centers = 2)$cluster)

#H3K27me3.clusters.wave1.genes <- unique(markers.mol[markers.mol$gene %in% rownames(heatmap.ls[[3]][[3]])[as.numeric(na.exclude(match(rownames(heatmap.ls[[3]][[3]]),names(a[a==1]))))],'closest_gene'])
#H3K27me3.clusters.wave2.genes <- unique(markers.mol[markers.mol$gene %in% rownames(heatmap.ls[[3]][[3]])[as.numeric(na.exclude(match(rownames(heatmap.ls[[3]][[3]]),names(a[a==2]))))],'closest_gene'])

H3K27me3.clusters.wave1.genes <- markers.mol[markers.mol$gene %in% rownames(heatmap.ls[[3]][[3]])[as.numeric(na.exclude(match(rownames(heatmap.ls[[3]][[3]]),names(a[a==1]))))],'closest_gene']
H3K27me3.clusters.wave2.genes <- markers.mol[markers.mol$gene %in% rownames(heatmap.ls[[3]][[3]])[as.numeric(na.exclude(match(rownames(heatmap.ls[[3]][[3]]),names(a[a==2]))))],'closest_gene']

# ens.gene <- genes(EnsDb.Mmusculus.v79)
# seqlevelsStyle(ens.gene) <- "UCSC"
# H3K27me3.clusters.wave1.genes <- unique(ens.gene[subjectHits(nearest(StringToGRanges(names(a[a==1])),ens.gene,select='all'))]$symbol)
# H3K27me3.clusters.wave2.genes <- unique(ens.gene[subjectHits(nearest(StringToGRanges(names(a[a==2])),ens.gene,select='all'))]$symbol)

#H3K27me3.clusters.wave1.genes <- gene_peak_table[match(gene_peak_table$query_region, names(a[a==1]),nomatch=0),'symbol']
#H3K27me3.clusters.wave2.genes <- gene_peak_table[match(gene_peak_table$query_region, names(a[a==2]),nomatch=0),'symbol']

# scRNA <- readRDS(file='/data/proj/GCB_MB/single-cell-CUT-Tag/nbiotech_paper/analysis/results/Sten_RNA/clustering/01.clustering_20000cells.Rds')
# scRNA <- scRNA[,scRNA$TaxonomyRank4 != 'Enteric neurons']
# scRNA <- scRNA[,scRNA$TaxonomyRank4 != 'Enteric glia']

scRNA.genes.scaled <- rownames(GetAssayData(object = scRNA,slot = 'scale.data',assay = 'RNA'))

r        <- match(H3K27me3.clusters.wave1.genes,scRNA.genes.scaled)
r        <- r[!is.na(r)]
r        <- scRNA.genes.scaled[r]
score.r  <- colMeans(GetAssayData(object = scRNA,slot = 'scale.data',assay = 'RNA')[r,])

scRNA <- AddMetaData(object = scRNA,metadata = score.r,col.name='Wave1.score')
p1 <- VlnPlot(scRNA,'Wave1.score', group.by='TaxonomyRank4',pt.size = 0) + NoLegend()
####

r        <- match(H3K27me3.clusters.wave2.genes,scRNA.genes.scaled)
r        <- r[!is.na(r)]
r        <- scRNA.genes.scaled[r]
score.r  <- colMeans(GetAssayData(object = scRNA,slot = 'scale.data',assay = 'RNA')[r,])

scRNA <- AddMetaData(object = scRNA,metadata = score.r,col.name='Wave2.score')
p2 <- VlnPlot(scRNA,'Wave2.score', group.by='TaxonomyRank4',pt.size = 0) + NoLegend()


# scRNA <- AddMetaData(object = scRNA,metadata = log2(scRNA$Wave2.score+0.0001/scRNA$Wave1.score+0.0001), col.name='Wave.ratio')
# p3 <- VlnPlot(scRNA,'Wave.ratio', group.by='TaxonomyRank4',pt.size=0) + NoLegend()
# p3

p1+p2 

p <- DimPlot(scRNA,group.by='Taxonomy_group',label=TRUE) + NoLegend()
p3 <- FeaturePlot(scRNA,'Wave1.score',max.cutoff=quantile(scRNA$Wave1.score,0.90)) + scale_color_viridis() + ggtitle('Wave2 genes')
p4 <- FeaturePlot(scRNA,'Wave2.score',max.cutoff=quantile(scRNA$Wave2.score,0.90)) + scale_color_viridis() + ggtitle('Wave1 genes')

p + p4 + p3

pdf('/data/proj/GCB_MB/temp.pdf',width=20,height=10)
grid.arrange(p,p4,p3)
dev.off()

#######################################################################################
# Which enteric glia genes contribute the most? 
r <- match(unique(H3K27me3.clusters.wave1.genes),rownames(scRNA))
r <- r[!is.na(r)]
barplot(sort(rowMeans(GetAssayData(scRNA)[r,which(scRNA$ClusterName == 'OPC')]),decreasing=TRUE))
head(sort(rowMeans(GetAssayData(scRNA)[r,which(scRNA$ClusterName == 'OPC')]),decreasing=TRUE),50)

#######################################################################################
r <- match(rownames(scRNA),unique(H3K27me3.clusters.wave2.genes))
r <- r[!is.na(r)]
barplot(sort(rowSums(GetAssayData(scRNA)[r,which(scRNA$ClusterName == 'COP1')]),decreasing=TRUE))
head(sort(rowSums(GetAssayData(scRNA)[r,which(scRNA$ClusterName == 'COP1')]),decreasing=TRUE),50)

################################# GO not working yet
library(goseq)
library(clusterProfiler)
library(org.Mm.eg.db)

ego1 <- enrichGO(gene          = H3K27me3.clusters.wave1.genes,
                 universe      = rownames(seurat.ls[[2]][['GA']]),
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = FALSE)

ego2 <- enrichGO(gene          = H3K27me3.clusters.wave2.genes,
                 universe      = rownames(seurat.ls[[2]][['GA']]),
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = FALSE)

ego3 <- gseDO(H3K27me3.clusters.wave1.genes)

#######################################################################################
marques <- readRDS(file='/data/proj/GCB_MB/single-cell-CUT-Tag/nbiotech_paper/analysis/results/marques_RNA/clustering/01.clustering.Rds')

r <- match(rownames(marques),H3K27me3.clusters.wave1.genes)
r <- r[!is.na(r)]

marques <- AddMetaData(object = marques,metadata = colMeans(GetAssayData(object = marques,slot = 'data',assay = 'RNA')[r,]),col.name='Wave1.score')
p1 <- VlnPlot(marques,'Wave1.score', group.by='cell_type',pt.size = 0,sort = 'decreasing') + NoLegend()

####

r <- match(rownames(marques),H3K27me3.clusters.wave2.genes)
r <- r[!is.na(r)]

marques <- AddMetaData(object = marques,metadata = colMeans(GetAssayData(object = marques,slot = 'data',assay = 'RNA')[r,]),col.name='Wave2.score')
p2 <- VlnPlot(marques,'Wave2.score', group.by='cell_type',pt.size=0,sort='decreasing') + NoLegend()


# scRNA <- AddMetaData(object = scRNA,metadata = log2(scRNA$Wave2.score+0.0001/scRNA$Wave1.score+0.0001), col.name='Wave.ratio')
# p3 <- VlnPlot(scRNA,'Wave.ratio', group.by='TaxonomyRank4',pt.size=0) + NoLegend()
# p3

p1 <- FeaturePlot(marques,'Wave1.score',max.cutoff=quantile(marques$Wave1.score,0.9)) + scale_color_viridis() + ggtitle('Wave2 genes')
p2 <- FeaturePlot(marques,'Wave2.score',max.cutoff=quantile(marques$Wave2.score,0.9)) + scale_color_viridis() + ggtitle('Wave1 genes')
p2+p1 

###################################





#######################################################################################
wave1.matrix.H3K27ac <- FeatureMatrix(fragments = Fragments(seurat.ls[[2]]),
                                      features = Extend(StringToGRanges(names(a[a==1])),upstream=10000,downstream=10000),
                                      cells = Cells(seurat.ls[[2]]))

wave2.matrix.H3K27ac <- FeatureMatrix(fragments = Fragments(seurat.ls[[2]]),
                                      features = Extend(StringToGRanges(names(a[a==2])),upstream=10000,downstream=10000),
                                      cells = Cells(seurat.ls[[2]]))


seurat.temp            <- CreateSeuratObject(counts = wave1.matrix.H3K27ac,assay = 'wave1')
seurat.temp[['wave2']] <- CreateAssayObject(counts=wave2.matrix.H3K27ac)

seurat.temp            <- SetIdent(object = seurat.temp,cells = names(seurat.ls[[2]]$idents_short),value = seurat.ls[[2]]$idents_short)

seurat.temp            <- AddMetaData(object = seurat.temp, colSums(wave1.matrix.H3K27ac),col.name = 'wave1_H3K27ac')
seurat.temp            <- AddMetaData(object = seurat.temp, colSums(wave2.matrix.H3K27ac),col.name = 'wave2_H3K27ac')

p1 <- VlnPlot(seurat.temp,features = 'wave1_H3K27ac',pt.size = 0.01,sort = 'decreasing')
p2 <- VlnPlot(seurat.temp,features = 'wave2_H3K27ac',pt.size = 0.01,sort = 'decreasing')

p1+p2

###################################################
df1 <-data.frame('score'    = colSums(heatmap.ls[['H3K27me3']][['ATAC']]), 
                 'pt'       = seurat.olg$pt[colnames(heatmap.ls[['ATAC']][['ATAC']])],
                 'modality' = 'ATAC')
p1 <- ggplot(data=df1,aes(x=pt,y=score)) + geom_point()

df2 <-data.frame('score'    = colSums(heatmap.ls[['H3K27me3']][['H3K27ac']]), 
                 'pt'       = seurat.olg$pt[colnames(heatmap.ls[['H3K27ac']][['H3K27ac']])],
                 'modality' = 'H3K27ac')
p2 <- ggplot(data=df2,aes(x=pt,y=score)) + geom_point()

df3 <-data.frame('score'    = colSums(heatmap.ls[['H3K27me3']][['H3K27me3']]), 
                 'pt'       = seurat.olg$pt[colnames(heatmap.ls[['H3K27me3']][['H3K27me3']])],
                 'modality' = 'H3K27me3')
p3 <- ggplot(data=df3,aes(x=pt,y=score)) + geom_point()

#grid.arrange(p1,p2,p3)

df.merged <- purrr::reduce(list(df1,df2,df3), rbind)
p <- ggplot(data=df.merged,aes(x=pt,y=score,col=modality)) + geom_point() + geom_smooth(method='loess')
p

#########################################################
df1 <-data.frame('score'    = colSums(heatmap.ls[['H3K27ac']][['ATAC']]), 
                 'pt'       = seurat.olg$pt[colnames(heatmap.ls[['H3K27ac']][['ATAC']])],
                 'modality' = 'ATAC')
p1 <- ggplot(data=df1,aes(x=pt,y=score)) + geom_point()

df2 <-data.frame('score'    = colSums(heatmap.ls[['H3K27ac']][['H3K27ac']]), 
                 'pt'       = seurat.olg$pt[colnames(heatmap.ls[['H3K27ac']][['H3K27ac']])],
                 'modality' = 'H3K27ac')
p2 <- ggplot(data=df2,aes(x=pt,y=score)) + geom_point()

df3 <-data.frame('score'    = colSums(heatmap.ls[['H3K27ac']][['H3K27me3']]), 
                 'pt'       = seurat.olg$pt[colnames(heatmap.ls[['H3K27me3']][['H3K27me3']])],
                 'modality' = 'H3K27me3')
p3 <- ggplot(data=df3,aes(x=pt,y=score)) + geom_point()

#grid.arrange(p1,p2,p3)

df.merged <- purrr::reduce(list(df1,df2,df3), rbind)
p <- ggplot(data=df.merged,aes(x=pt,y=score,col=modality)) + geom_point() + geom_smooth(method='loess')
p




########################################################################

df1 <-data.frame('score'    = colSums(heatmap.ls[['ATAC']][['ATAC']]), 
                 'pt'       = seurat.olg$pt[colnames(heatmap.ls[['ATAC']][['ATAC']])],
                 'modality' = 'ATAC')
p1 <- ggplot(data=df1,aes(x=pt,y=score)) + geom_point()

df2 <-data.frame('score'    = colSums(heatmap.ls[['ATAC']][['H3K27ac']]), 
                 'pt'       = seurat.olg$pt[colnames(heatmap.ls[['ATAC']][['H3K27ac']])],
                 'modality' = 'H3K27ac')
p2 <- ggplot(data=df2,aes(x=pt,y=score)) + geom_point()

df3 <-data.frame('score'    = colSums(heatmap.ls[['ATAC']][['H3K27me3']]), 
                 'pt'       = seurat.olg$pt[colnames(heatmap.ls[['ATAC']][['H3K27me3']])],
                 'modality' = 'H3K27me3')
p3 <- ggplot(data=df3,aes(x=pt,y=score)) + geom_point()

#grid.arrange(p1,p2,p3)

df.merged <- purrr::reduce(list(df1,df2,df3), rbind)
p <- ggplot(data=df.merged,aes(x=pt,y=score,col=modality)) + geom_point() + geom_smooth(method='loess')
p
########################################################################




colsums.heatmap.ls <- list()
for(i in modalities){
  colsums.heatmap.ls[[i]] <- list()
  for(x in modalities){
    colsums.heatmap.ls[[i]][[x]] <- data.frame('score' = colSums(heatmap.ls[[i]][[x]]),
                                               'pt' = seurat.olg$pt[colnames(heatmap.ls[[i]][[x]])],
                                               'markers' = i,
                                               'signal' = x)
  }
}



p1 <- ggplot(data=colsums.heatmap.ls[['ATAC']][['ATAC']], aes(x=pt,y=score)) + geom_point() + theme_minimal() + geom_smooth(method = 'loess')
p2 <- ggplot(data=colsums.heatmap.ls[['ATAC']][['H3K27ac']], aes(x=pt,y=score)) + geom_point() + theme_minimal()  + geom_smooth(method = 'loess')
p3 <- ggplot(data=colsums.heatmap.ls[['ATAC']][['H3K27me3']], aes(x=pt,y=score)) + geom_point() + theme_minimal()  + geom_smooth(method = 'loess')
p1 + p2 + p3



p1 <- ggplot(data=colsums.heatmap.ls[['H3K27ac']][['ATAC']], aes(x=pt,y=score)) + geom_point() + theme_minimal() + geom_smooth(method = 'loess')
p2 <- ggplot(data=colsums.heatmap.ls[['H3K27ac']][['H3K27ac']], aes(x=pt,y=score)) + geom_point() + theme_minimal()  + geom_smooth(method = 'loess')
p3 <- ggplot(data=colsums.heatmap.ls[['H3K27ac']][['H3K27me3']], aes(x=pt,y=score)) + geom_point() + theme_minimal()  + geom_smooth(method = 'loess')
grid.arrange(p1,p2,p3,ncol=1)


p1 <- ggplot(data=colsums.heatmap.ls[['H3K27me3']][['ATAC']], aes(x=pt,y=score)) + geom_point() + theme_minimal() + geom_smooth(method = 'loess')
p2 <- ggplot(data=colsums.heatmap.ls[['H3K27me3']][['H3K27ac']], aes(x=pt,y=score)) + geom_point() + theme_minimal()  + geom_smooth(method = 'loess')
p3 <- ggplot(data=colsums.heatmap.ls[['H3K27me3']][['H3K27me3']], aes(x=pt,y=score)) + geom_point() + theme_minimal()  + geom_smooth(method = 'loess')
p1 + p2 + p3


colsums.ls <- lapply(colsums.heatmap.ls,function(x){purrr::reduce(x,rbind)})
p1 <- ggplot(data=colsums.ls[[1]], aes(x=pt,y=score,col=signal)) + geom_point() + geom_smooth(method = 'loess')
p2 <- ggplot(data=colsums.ls[[2]], aes(x=pt,y=score,col=signal)) + geom_point() + geom_smooth(method = 'loess')
p3 <- ggplot(data=colsums.ls[[3]], aes(x=pt,y=score,col=signal)) + geom_point() + geom_smooth(method = 'loess')

p1
p2
p3


grid.arrange(p1,p2,p3)


seurat.ls <- readRDS(file='/data/proj/GCB_MB/bcd_CT/single-cell/results/multiple_modalities/ATAC_H3K27ac_H3K27me3/seurat_multimodal/peaks/Seurat_object.Rds')
seurat.oligo.ls <- lapply(seurat.ls[1:3],function(x){
  x.oligo <- x[,x$idents_short %in% c('OPC','MOL')]
  x.oligo <- AddMetaData(object=x.oligo,metadata = seurat.olg$pt,col.name = 'pt')
})

# Cicero
run_cicero_MB <- function(seurat){

  cds        <- as.CellDataSet(seurat)
  cds.cicero <- make_cicero_cds(cds, reduced_coordinates = reducedDimS(cds))

  genome <- seqlengths(seurat)
  genome.df <- data.frame("chr" = names(genome), "length" = genome)

  conns <- run_cicero(cds.cicero, genomic_coords = genome.df, sample_num = 100)
  ccans <- generate_ccans(conns)
  return(ccans)
}

ccans.ls <- lapply(seurat.oligo.ls,run_cicero_MB)


# Get common cells
cells <- lapply(seurat.oligo.ls,Cells)
cells <- purrr::reduce(cells,intersect)

cicero.regions.matrix.ls <- lapply(ccans.ls, function(ccan){
  regions <- StringToGRanges(ccan$Peak)
  cc_score.df <- data.frame('cell' = cells,
                            'pt' = seurat.olg$pt[cells],
                            row.names = cells)
  feature_matrix_ls <- list()
  r.ls <- lapply(seurat.oligo.ls,function(s){
    feature_matrix <- Signac::FeatureMatrix(fragments = Fragments(s),
                                            features = regions,
                                            cells = cells)
  })
})

cicero.regions.scores.ls <- list()
for(i in names(cicero.regions.matrix.ls)){
  cicero.regions.scores.ls[[i]] <- list()
  print(i)
  for(j in names(cicero.regions.matrix.ls[[i]])){
    print(j)
    cicero.regions.scores.ls[[i]][[j]] <- data.frame('cell' = cells,
                                                 'pt' = seurat.olg$pt[cells],
                                                  row.names = cells)
    for(cc in unique(ccans.ls[[i]]$CCAN)){
      regions  <- ccans.ls[[i]][ccans.ls[[i]]$CCAN == cc ,'Peak']
      cc.score <- cicero.regions.matrix.ls[[i]][[j]][regions,]
      cc.score <- colSums(cc.score)
      cicero.regions.scores.ls[[i]][[j]][,paste0('CC_',cc)] <- cc.score
    }
  }
}

pt <- sort(seurat.olg[,cells]$pt)

CC.scores.df <- lapply(cicero.regions.scores.ls[[2]],function(scores){
  df <- scores
  df <- df[!is.na(df$pt),]
  df.meta <- df[,1:2]
  df <- as.matrix(df[,-1:-2])

  df <- scale(df)
  df <- cbind(df.meta,df)

  correlation <- list()
  for(cc in colnames(df)[3:dim(df)[2]]){correlation[[cc]] <- lm(df[,cc]~ df$pt)$coefficients[2] }
  correlation        <- unlist(correlation,recursive = TRUE)
  names(correlation) <- unlist(lapply(strsplit(names(correlation),'\\.'),'[',1))

  print(head(sort(correlation,decreasing=TRUE),20))
  return(df)
})
# CC_1021   CC_1670   CC_2346   CC_2101   CC_2019   CC_2106   CC_1555   CC_1663   CC_1913    CC_759   CC_2238   CC_1814   CC_1007 
# 0.3058981 0.2724308 0.2661198 0.2450364 0.2446587 0.2429624 0.2403327 0.2402232 0.2385427 0.2379687 0.2347453 0.2330394 0.2325070 
# CC_893   CC_1443   CC_2008    CC_969   CC_1683   CC_1537    CC_933 
# 0.2320496 0.2313055 0.2285762 0.2261408 0.2243296 0.2238814 0.2236666 

to.plot.ls <- lapply(names(CC.scores.df),function(x,CC){
  r <- data.frame('pt' = CC.scores.df[[x]][,'pt'],
                  'score' = CC.scores.df[[x]][,'CC_2106'],
                  'modality'=x)
  r
})
to.plot <- purrr::reduce(to.plot.ls,rbind)
ggplot(data=to.plot, aes(x=pt,y=score,col=modality)) + geom_point() + geom_smooth(method='loess')

p.ls <- lapply(names(head(sort(correlation,decreasing=TRUE),200)), function(x){
  r <- df[names(pt),x]
})

to.plot <- purrr::reduce(p.ls,rbind)
head(to.plot.ls[[1]][1:10,1:10])

plot(colSums(to.plot.ls[[2]]))
pheatmap(to.plot,cluster_cols = F,cluster_rows = F)#,color = c(rep('black',50),viridis_pal()(50),rep('white',50)),scale='row')


############################################### scRNA-seq






