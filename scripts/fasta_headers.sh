#!/bin/bash


analysis=analysis_name

mkdir $analysis
mkdir $analysis/fasta_out
mkdir $analysis/fasta_out/selected_ogs
mkdir $analysis/fasta_out/selected_ogs/original_headers
mkdir $analysis/fasta_out/selected_ogs/sp_only_headers


filelist=$analysis\_ogs.txt
data=/path/to/filtered/aligned/backtranslated/files/$analysis\_out/filtered_out/

ogs=$(cat $filelist)

for f in $ogs;
do
	file=$(find $data -name "aligned_mafft_${f}_filtered.fa")
	cp $file $analysis/fasta_out/selected_ogs/original_headers
	sed 's/^>.*_/>/; s/_[0-9]+ *$//' $analysis/fasta_out/selected_ogs/original_headers/aligned_mafft_"$f"_filtered.fa > $analysis/fasta_out/selected_ogs/sp_only_headers/aligned_mafft_"$f"_filtered_sp_only.fa


done;

