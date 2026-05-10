#!/bin/bash


ANALYSIS=analysis_name

INFOLDER=$(ls -d /path/to/backtranslated/sequences/$ANALYSIS\_out/backtranslated_out/hits_per_*_species)

rm -r $ANALYSIS\_out
mkdir $ANALYSIS\_out

OUTFOLDER=/path/to/output/$ANALYSIS\_out


for i in $INFOLDER;
do


	OGS=$(ls "$i"/*.fa)

	for o in $OGS;
	do
		n=$(basename "$o" .fa_aligned_mafft_pal2nal_out.fa)
		seqkit grep -r -p .*Prang.* "$o" > $OUTFOLDER/"$n"_Prang_OGs.fa
		sed "s/>.*/&-${n}/" $OUTFOLDER/"$n"_Prang_OGs.fa > $OUTFOLDER/"$n"_Prang_OGs_2.fa
		cat $OUTFOLDER/*_Prang_OGs_2.fa > /path/to/output/$ANALYSIS\_Prang_OGs.fa
	done

done
