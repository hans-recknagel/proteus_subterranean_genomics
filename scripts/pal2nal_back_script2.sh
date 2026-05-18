#!/bin/bash

#Input is the path to folder containing all filtered orthogroups and their mafft alignments.
#The rest of the paths should be set accordingly so the script can find the AA alignments and compare them to the nucleotide sequences. 

INFOLDER=/path/to/analysis_name

DIRS=$(ls -d $INFOLDER | cat)

cd $INFOLDER

NSP=$(ls -d filtered_out/*species | sed -e s/[^0-9]//g)
mkdir pal2nal_out/

for n in $NSP;
do
	mafftdir=~/path/to/$INFOLDER/backtranslated_out/hits_per_"$n"_species_aligned_mafft/
	muscledir=~/path/to/$INFOLDER/backtranslated_out/hits_per_"$n"_species_aligned_muscle_stable/
	if [ -d "$mafftdir" ];
	then
		mkdir pal2nal_out/hits_per_"$n"_species_mafft_pal2nal
	fi
	if [ -d "$muscledir" ];
	then
		mkdir pal2nal_out/hits_per_"$n"_species_muscle_stable_pal2nal
	fi
	cd backtranslated_out/
	cd hits_per_"$n"_species/
	FILES=$(ls)
	for f in $FILES;
	do
		filename="${f%.*}"
		if [ -d "$mafftdir" ];
		then
			/path/to/pal2nal.pl \
				/path/to/$INFOLDER/filtered_out/hits_per_"$n"_species_aligned_mafft/aligned_mafft_"$filename"_filtered.fa \
				"$f" \
				-output fasta > ../../pal2nal_out/hits_per_"$n"_species_mafft_pal2nal/"$f"_aligned_mafft_pal2nal_out.fa
		fi
		if [ -d "$muscledir" ];
		then
			/path/to/pal2nal.pl \
				/path/to/$INFOLDER/filtered_out/hits_per_"$n"_species_aligned_muscle_stable/aligned_muscle_stable_"$filename"_filtered.fa \
				"$f" \
				-output fasta > ../../pal2nal_out/hits_per_"$n"_species_muscle_stable_pal2nal/"$f"_aligned_muscle_stable_pal2nal_out.fa

		fi

	done;
	cd ../..
done;
