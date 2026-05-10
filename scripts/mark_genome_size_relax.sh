#!/bin/bash

files=$(ls best_trees_relax/genome_size/*nwk)
#cat $files

for f in $files;
do
	sed -i 's/Plwal/Plwal{Foreground}/g' "$f"
	sed -i 's/Ammex/Ammex{Foreground}/g' "$f"
	sed -i 's/Prang/Prang{Foreground}/g' "$f"
	sed -i 's/Bubuf/Bubuf{Foreground}/g' "$f"
	sed -i 's/Nefor/Nefor{Foreground}/g' "$f"
	sed -i 's/Prann/Prann{Foreground}/g' "$f"
	sed -i 's/Sctor/Sctor{Foreground}/g' "$f"
	sed -i 's/Asmes/Asmes{Background}/g' "$f"
	sed -i 's/Asmec/Asmec{Background}/g' "$f"
	sed -i 's/Darer/Darer{Background}/g' "$f"
	sed -i 's/Hosap/Hosap{Background}/g' "$f"
	sed -i 's/Leocu/Leocu{Background}/g' "$f"
	sed -i 's/Tagut/Tagut{Background}/g' "$f"
	sed -i 's/Ornil/Ornil{Background}/g' "$f"
	sed -i 's/Camil/Camil{Background}/g' "$f"
	sed -i 's/Enpus/Enpus{Background}/g' "$f"
	sed -i 's/Lacha/Lacha{Background}/g' "$f"
	sed -i 's/Menas/Menas{Background}/g' "$f"
	sed -i 's/Pemar/Pemar{Background}/g' "$f"
	sed -i 's/Pipip/Pipip{Background}/g' "$f"
	sed -i 's/Xetro/Xetro{Background}/g' "$f"
	sed -i 's/Zoviv/Zoviv{Background}/g' "$f"
done
