#!/bin/bash

files=$(ls best_trees_relax/longevity/*nwk)
#cat $files

for f in $files;
do
        sed -i 's/Acbae/Acbae{Foreground}/g' "$f"
        sed -i 's/Cacar/Cacar{Foreground}/g' "$f"
        sed -i 's/Prang/Prang{Foreground}/g' "$f"
        sed -i 's/Elmax/Elmax{Foreground}/g' "$f"
        sed -i 's/Hosap/Hosap{Foreground}/g' "$f"
        sed -i 's/Lacha/Lacha{Foreground}/g' "$f"
        sed -i 's/Myluc/Myluc{Foreground}/g' "$f"
        sed -i 's/Hegla/Hegla{Foreground}/g' "$f"
        sed -i 's/Ororc/Ororc{Foreground}/g' "$f"
        sed -i 's/Rhtyp/Rhtyp{Foreground}/g' "$f"
        sed -i 's/Tecar/Tecar{Foreground}/g' "$f"
	sed -i 's/Ancor/Ancor{Background}/g' "$f"
	sed -i 's/Botau/Botau{Background}/g' "$f"
	sed -i 's/Capor/Capor{Background}/g' "$f"
	sed -i 's/Darer/Darer{Background}/g' "$f"
	sed -i 's/Chsab/Chsab{Background}/g' "$f"
	sed -i 's/Chpic/Chpic{Background}/g' "$f"
	sed -i 's/Cocri/Cocri{Background}/g' "$f"
	sed -i 's/Ercal/Ercal{Background}/g' "$f"
	sed -i 's/Modom/Modom{Background}/g' "$f"
	sed -i 's/Ovari/Ovari{Background}/g' "$f"
	sed -i 's/Ranor/Ranor{Background}/g' "$f"
	sed -i 's/Sasal/Sasal{Background}/g' "$f"
	sed -i 's/Scvul/Scvul{Background}/g' "$f"
	sed -i 's/Tagut/Tagut{Background}/g' "$f"
done
