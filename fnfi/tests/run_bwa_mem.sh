#!/bin/bash

bwa_ref=$1  # e.g. ~/Documents/Data/db/hg38/hg38.fa
sample=$2   # e.g. ./sample.fastq
threads=$3

bwa mem -t$3 -P -a -p $bwa_ref $sample