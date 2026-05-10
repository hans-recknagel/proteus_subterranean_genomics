#!/bin/bash

analysis=analysis_name
hits=hits_per_(insert_number)_species

fastadir=/path/to/$analysis\_out/pal2nal_out/$hits\_mafft_pal2nal
treedir=/path/to/best_trees_absrel/$analysis
outdir=/path/to/absrel/$analysis

mkdir $outdir
mkdir $outdir/$hits


files=$(ls $fastadir/*fa)
for f in $files;
do
	nameonly=$(basename "$f")
	hyphy CPU=10 absrel \
		--alignment "$f" \
		--tree $treedir/RAxML_bestTree."$nameonly".nwk \
		--output $outdir/$hits/"$nameonly".ABSREL.json
done




