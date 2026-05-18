#!bin/bash

#Inputs are analysis name, number of species, pal2nal output alignments, RAxML trees, and the path to relax results.


analysis=anlysis_name


treedir=/path/to/best_trees_relax/$analysis



files=$(find relax/$analysis* -name '*.json' -size 0)


for f in $files;
do
        nameonly=$(basename "$f" | sed 's/\.0[1-5].*//')
        hits1=$(dirname "$f")
        hits=$(basename $hits1)

        fastadir=/data/home/lukam/proteus_genome/orthologs/backtranslatingOFresults/$analysis\_out/pal2nal_out/$hits\_mafft_pal2nal

        hyphy CPU=10 ENV="TOLERATE_NUMERICAL_ERRORS=1;" relax \
                --alignment $fastadir/$nameonly \
                --tree $treedir/RAxML_bestTree."$nameonly".nwk \
                --test "Foreground" \
                --reference "Background" \
                --starting-points 10 \
                --grid-size 1000 \
                --output "$f"

done





