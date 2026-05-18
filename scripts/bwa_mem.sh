#!/bin/bash

#The input is the split genome and results from the filter_for_prang.sh script.

REF=/path/to/genome/genome.fasta

FASTA=$(ls *_Prang_OGs.fa)
OUTDIR=bwa_mem_out

rm -r $OUTDIR
mkdir $OUTDIR


for f in $FASTA;
do
	name=$(basename "$f" _Prang_OGs.fa)
	bwa mem -t 10 $REF "$f" > $OUTDIR/$name.bwa.aln.sam

done
