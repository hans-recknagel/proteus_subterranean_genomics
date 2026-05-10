#!/bin/bash

INFOLDER=infolder_name
INFOLDER2=/absolute/path/to/$INFOLDER

TXT=$(ls $INFOLDER2)

for t in $TXT;
do
	if [[ "$t" == *"cave"* ]];
	then
		ANALISYS=cave_species
	elif [[ "$t" == *"eyes"* ]]
	then
		ANALISYS=eye_loss
        elif [[ "$t" == *"geno"* ]]
        then
		ANALISYS=genome_size
        elif [[ "$t" == *"long"* ]]
        then
		ANALISYS=longevity
        elif [[ "$t" == *"pigm"* ]]
        then
		ANALYSIS=pigmentation_loss
	fi;

	OGS=$(cat $INFOLDER/"$t" | tail -n +2 | awk '{print $1}')

	for o in $OGS;
	do
		OG=$(find /path/to/orthofinder_results/$ANALISYS/Orthogroup_Sequences/ -name "*${o}*" -print)

		seqkit grep -r -p .*Prang.* $OG >> "$t"_Prang_OGs.fa

		sed "s/>.*/&-${o}/" "$t"_Prang_OGs.fa > "$t"_Prang_OGs_2.fa

	done

done
