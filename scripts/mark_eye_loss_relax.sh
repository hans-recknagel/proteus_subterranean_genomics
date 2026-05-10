#!/bin/bash

files=$(ls best_trees_relax/eye_loss/*nwk)
#cat $files

for f in $files;
do

        sed -i 's/Asmec/Asmec{Foreground}/g' "$f"
        sed -i 's/Chasi/Chasi{Foreground}/g' "$f"
        sed -i 's/Luden/Luden{Foreground}/g' "$f"
        sed -i 's/Cocri/Cocri{Foreground}/g' "$f"
        sed -i 's/Geser/Geser{Foreground}/g' "$f"
        sed -i 's/Hegla/Hegla{Foreground}/g' "$f"
        sed -i 's/Miuni/Miuni{Foreground}/g' "$f"
        sed -i 's/Nagal/Nagal{Foreground}/g' "$f"
        sed -i 's/Prang/Prang{Foreground}/g' "$f"
        sed -i 's/Sians/Sians{Foreground}/g' "$f"
        sed -i 's/Trros/Trros{Foreground}/g' "$f"
	sed -i 's/Ammex/Ammex{Background}/g' "$f"
	sed -i 's/Ammel/Ammel{Background}/g' "$f"
	sed -i 's/Asmes/Asmes{Background}/g' "$f"
	sed -i 's/Darer/Darer{Background}/g' "$f"
	sed -i 's/Gamor/Gamor{Background}/g' "$f"
	sed -i 's/Mymur/Mymur{Background}/g' "$f"
	sed -i 's/Plwal/Plwal{Background}/g' "$f"
	sed -i 's/Sigra/Sigra{Background}/g' "$f"
	sed -i 's/Thmac/Thmac{Background}/g' "$f"
	sed -i 's/Trtib/Trtib{Background}/g' "$f"
	sed -i 's/Bubuf/Bubuf{Background}/g' "$f"
	sed -i 's/Chlan/Chlan{Background}/g' "$f"
	sed -i 's/Ectel/Ectel{Background}/g' "$f"
	sed -i 's/Elmax/Elmax{Background}/g' "$f"
	sed -i 's/Ereur/Ereur{Background}/g' "$f"
	sed -i 's/Jajac/Jajac{Background}/g' "$f"
	sed -i 's/Ornil/Ornil{Background}/g' "$f"
	sed -i 's/Orcun/Orcun{Background}/g' "$f"
	sed -i 's/Ranor/Ranor{Background}/g' "$f"
	sed -i 's/Soara/Soara{Background}/g' "$f"
	sed -i 's/Sirhi/Sirhi{Background}/g' "$f"
done
