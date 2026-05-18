#!/bin/bash

#Input is the path to the fasta files from the fasta_headers_BF_*.sh scripts, name of the analysis, and an output path.

analysis=analysis_name

fastapath=$analysis/fasta_out/selected_ogs/sp_only_headers/

files=$(ls $fastapath/*fa)
outpath=$analysis/phylip
mkdir $outpath

for f in $files;
do
	name=$(basename "$f" .fa)
	seqret -sformat1 pearson -osformat2 phylipnon -sequence "$f" -outseq $outpath/$name.phylip
done;

