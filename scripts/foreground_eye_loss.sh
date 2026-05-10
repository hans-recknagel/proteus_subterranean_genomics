#!/bin/bash

files=$(ls best_trees_absrel/eye_loss/*nwk)
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
done
