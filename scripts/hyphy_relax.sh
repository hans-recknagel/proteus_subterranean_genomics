#!/bin/bash

analysis=analysis_name
hits=hits_per_(insert number of species)_species
version=01

fastadir=/path/to/$analysis\_out/pal2nal_out/$hits\_mafft_pal2nal
treedir=/path/to/best_trees_relax/$analysis
outdir=/path/to/relax/$analysis.$version.GITHUBFIXES

mkdir $outdir
mkdir $outdir/$hits


files=$(ls $fastadir/*fa)


for f in $files;
do
        nameonly=$(basename "$f")
        hyphy CPU=10 ENV="TOLERATE_NUMERICAL_ERRORS=1;" relax \
                --alignment "$f" \
                --tree $treedir/RAxML_bestTree."$nameonly".nwk \
                --test "Foreground" \
                --reference "Background" \
                --starting-points 10 \
                --grid-size 1000 \
                --output $outdir/$hits/"$nameonly".$version.GITHUBFIXES.RELAX.json

done

