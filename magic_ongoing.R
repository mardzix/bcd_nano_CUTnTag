---
title: "R Notebook"
output: html_notebook
---

```{r}
library(Rmagic)
library(Seurat)
library(Signac)
set.seed(1234)
```

```{r}
seurat.ls <- readRDS(file='/data/proj/GCB_MB/bcd_CT/single-cell/results/multiple_modalities/ATAC_H3K27ac_H3K27me3/seurat_multimodal/peaks/Seurat_object.Rds')
```

```{r}
# imputed <- magic(seurat.ls[[1]],genes = 'chr9-69754058-69772974')
# FeaturePlot(imputed,'chr9-69754058-69772974') 
# 
# seurat.magic.ls <- lapply(seurat.ls,magic)
# lapply(seurat.magic.ls,function(x){DefaultAssay(x) <- 'MAGIC_peaks'})
# 
# FeaturePlot(seurat.magic.ls[[1]],'chr9-69754058-69772974') 

peaks.common <- rtracklayer::import(con = '/data/proj/GCB_MB/bcd_CT/single-cell/results/benchmarks/peaks/specificity_benchmark/peaks_merged.bed')

##############
peaks.ls <- lapply(seurat.ls[1:3],function(x){StringToGRanges(rownames(x))})
peaks.merged <- purrr::reduce(peaks.ls,c)
peaks.merged <- GenomicRanges::reduce(peaks.merged)
##############


matrix.new.ls <- lapply(seurat.ls[1:3],function(x){
  FeatureMatrix(fragments = Fragments(x),features = peaks.merged,cells = Cells(x))
})

seurat.common.ls <- lapply(names(seurat.ls)[1:3], function(x){
  CreateSeuratObject(counts = matrix.new.ls[[x]],assay = x)
})

seurat.common.magic.ls <- lapply(seurat.common.ls, magic)

```

