#!/bin/bash

#Inputs are the analysis name, list of orthogroups for that analysis and part to the tree folder.

analysis=analysis_name
list=$analysis\_ogs_all_OGonly.txt
ogs=$(cat $list)

treedir=/mnt/disk1/zoo/luka/proteus_genome/trees_input/$analysis

outdir=$analysis/originalTrees
mkdir $outdir

for o in $ogs;
do
	cp $treedir/RAxML_bestTree."$o".fa_aligned_mafft_pal2nal_out.fa.nwk $outdir/RAxML_bestTree."$o".fa_aligned_mafft_pal2nal_out.fa.nwk
done;
