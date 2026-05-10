#!/bin/bash

files=$(ls best_trees_absrel/genome_size/*nwk)
#cat $files

for f in $files;
do
	sed -i 's/Plwal/Plwal{Foreground}/g' "$f"
	sed -i 's/Ammex/Ammex{Foreground}/g' "$f"
	sed -i 's/Prang/Prang{Foreground}/g' "$f"
	sed -i 's/Prann/Prann{Foreground}/g' "$f"
	sed -i 's/Sctor/Sctor{Foreground}/g' "$f"
	sed -i 's/Bubuf/Bubuf{Foreground}/g' "$f"
	sed -i 's/Nefor/Nefor{Foreground}/g' "$f"
done
