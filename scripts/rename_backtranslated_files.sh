#!/bin/bash

#Inputs are backtranslated files from backtranslation_cleaned_files_exp_con_results_for_bwa.R 

INFOLDER=infolder_name
OUTFOLDER=infolder_name_new_headers

rm -r $OUTFOLDER
mkdir $OUTFOLDER

FOLDERS=$(ls -d $INFOLDER/*)

for f in $FOLDERS;
do

	foldername=$(basename "$f")
	mkdir $OUTFOLDER/$foldername
	OGS=$(ls "$f")
	for o in $OGS;
	do
		filename=$(basename "$o" .fa)
		sed "s/>/>${foldername}-${filename}-/" "$f"/"$o" > $OUTFOLDER/$foldername/$filename\_renamed.fa

	done

	cat $OUTFOLDER/$foldername/*_renamed.fa > $OUTFOLDER/$foldername\_all.fa

done
