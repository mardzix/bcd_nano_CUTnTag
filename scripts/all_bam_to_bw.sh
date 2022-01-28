#!/usr/bin/env bash

# $1 = path to folder with bam files [input]
# $2 = path to folder with bw files [output]
# $3 = nthreads

mkdir $2

ls $1"/"*.bam | while read line; do
  OUT_BW=$2"/"`basename ${line/.bam/.bw}`
  echo $line
  echo $OUT_BW
  samtools index $line
  bamCoverage -b $line -o $OUT_BW -p $3 --minMappingQuality 5 --binSize 50 --centerReads --smoothLength 250 --normalizeUsing RPKM --ignoreDuplicates --extendReads
  done

