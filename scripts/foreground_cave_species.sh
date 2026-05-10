#!/bin/bash

files=$(ls cave_species/*nwk)
#cat $files

for f in $files;
do
	sed -i 's/Asmec/Asmec{Foreground}/g' "$f"
	sed -i 's/Luden/Luden{Foreground}/g' "$f"
	sed -i 's/Prang/Prang{Foreground}/g' "$f"
	sed -i 's/Sians/Sians{Foreground}/g' "$f"
	sed -i 's/Sirhi/Sirhi{Foreground}/g' "$f"
	sed -i 's/Trros/Trros{Foreground}/g' "$f"
done
