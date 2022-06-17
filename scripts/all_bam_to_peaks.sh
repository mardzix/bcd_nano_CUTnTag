#!/usr/bin/env bash

# $1 = path to folder with bam files [input]
# $2 = path to folder with peaks [output]

mkdir $2

ls $1"/"*.bam | while read line; do
  OUT=$2"/"`basename ${line/.bam/}`"/"
  NAME=`basename ${line/.bam/}`
  echo $line
  echo $OUT
  echo $NAME
  macs2 callpeak -t $line -g mm -f BAMPE -n $NAME --outdir $OUT --keep-dup=1 --llocal 100000 --cutoff-analysis --min-length 1000 --max-gap 1000
  done

