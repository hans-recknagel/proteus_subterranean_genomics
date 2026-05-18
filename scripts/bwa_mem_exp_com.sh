#!/bin/bash

#Inputs are the path to the renamed backtranslated results from script rename_backtranslated_files.sh, output path, and the split genome.

FOLDER=infolder_name_new_headers
OUTFOLDER=outfolder_name

REF=/path/to/ref/genome.fasta

rm -r $OUTFOLDER
mkdir $OUTFOLDER

FASTA=$(ls $FOLDER/*_all.fa)

for f in $FASTA;
do
	name=$(basename "$f" _all.fa)
	bwa mem -t 10 $REF "$f" > $OUTFOLDER/$name.bwa.aln.sam

done
