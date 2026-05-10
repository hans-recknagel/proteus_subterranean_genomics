#!/bin/bash

files=$(ls best_trees_relax/pigmentation_loss/*nwk)
#cat $files

for f in $files;
do
	sed -i 's/Ammex/Ammex{Foreground}/g' "$f"
	sed -i 's/Asmec/Asmec{Foreground}/g' "$f"
	sed -i 's/Hegla/Hegla{Foreground}/g' "$f"
	sed -i 's/Prang/Prang{Foreground}/g' "$f"
	sed -i 's/Luden/Luden{Foreground}/g' "$f"
	sed -i 's/Sirhi/Sirhi{Foreground}/g' "$f"
	sed -i 's/Sians/Sians{Foreground}/g' "$f"
	sed -i 's/Trros/Trros{Foreground}/g' "$f"
	sed -i 's/Ammel/Ammel{Background}/g' "$f"
	sed -i 's/Asmes/Asmes{Background}/g' "$f"
	sed -i 's/Darer/Darer{Background}/g' "$f"
	sed -i 's/Gamor/Gamor{Background}/g' "$f"
	sed -i 's/Hosap/Hosap{Background}/g' "$f"
	sed -i 's/Leocu/Leocu{Background}/g' "$f"
	sed -i 's/Mymur/Mymur{Background}/g' "$f"
	sed -i 's/Plwal/Plwal{Background}/g' "$f"
	sed -i 's/Sigra/Sigra{Background}/g' "$f"
	sed -i 's/Tagut/Tagut{Background}/g' "$f"
	sed -i 's/Thmac/Thmac{Background}/g' "$f"
	sed -i 's/Xetro/Xetro{Background}/g' "$f"
	sed -i 's/Zoviv/Zoviv{Background}/g' "$f"
	sed -i 's/Trtib/Trtib{Background}/g' "$f"
	sed -i 's/Chlan/Chlan{Background}/g' "$f"
	sed -i 's/Ornil/Ornil{Background}/g' "$f"
	sed -i 's/Scvul/Scvul{Background}/g' "$f"
done
