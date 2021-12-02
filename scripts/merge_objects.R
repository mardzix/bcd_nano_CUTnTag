#!/usr/bin/env Rscript

library(Seurat)
library(Signac)
library(argparse)
library(purrr)
set.seed(1234)

########### Arguments parser

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default='foo', 
                    help="seurat objects to be merged",nargs="+")
parser$add_argument("-o", "--output", type="character", default='foo', 
                    help="output file")
args <- parser$parse_args()
######################## End arguments parser

seurat.ls <- lapply(args$input,function(x){
  readRDS(file=x)
})

seurat.merged <- purrr::reduce(seurat.ls,merge)
saveRDS(object = seurat.merged,file = args$output)

