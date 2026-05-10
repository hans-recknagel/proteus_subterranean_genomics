#!/bin/bash

files=$(ls best_trees_absrel/longevity/*nwk)
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
done
