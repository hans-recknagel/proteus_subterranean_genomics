#!/bin/bash

#Input is the name of the analysis, and path to the filtered mafft aligned files.

analysis=analysis_name

mkdir $analysis
mkdir $analysis/fasta_out
mkdir $analysis/fasta_out/selected_ogs
mkdir $analysis/fasta_out/selected_ogs/original_headers
mkdir $analysis/fasta_out/selected_ogs/sp_only_headers

find  /path/to/filtered/aligned/backtranslated/files/$analysis\_out/filtered_out/hits_per_*_species_aligned_mafft -name 'aligned_mafft_*_filtered.fa' -type f > $analysis\_ogs_all.txt
grep -o 'OG[0-9][0-9][0-9][0-9][0-9][0-9][0-9]' $analysis\_ogs_all.txt > $analysis\_ogs_all_OGonly.txt

ogs=$(cat $analysis\_ogs_all.txt)

for f in $ogs;
do
	name=$(basename "$f" _filtered.fa)
	name2=$(basename "$f")
	cp "$f" $analysis/fasta_out/selected_ogs/original_headers
	sed 's/^>.*_/>/; s/_[0-9]+ *$//' $analysis/fasta_out/selected_ogs/original_headers/$name2 > $analysis/fasta_out/selected_ogs/sp_only_headers/$name\_filtered_sp_only.fa


done;

