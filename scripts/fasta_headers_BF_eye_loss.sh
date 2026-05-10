#!/bin/bash

files=$(ls eye_loss/fasta_out/selected_ogs/sp_only_headers/*fa)
#cat $files

for f in $files;
do
	sed -i 's/Asmec/Fo_Asmec/g' "$f"
        sed -i 's/Chasi/Fo_Chasi/g' "$f"
        sed -i 's/Luden/Fo_Luden/g' "$f"
        sed -i 's/Cocri/Fo_Cocri/g' "$f"
        sed -i 's/Geser/Fo_Geser/g' "$f"
        sed -i 's/Hegla/Fo_Hegla/g' "$f"
        sed -i 's/Miuni/Fo_Miuni/g' "$f"
        sed -i 's/Nagal/Fo_Nagal/g' "$f"
        sed -i 's/Prang/Fo_Prang/g' "$f"
        sed -i 's/Sians/Fo_Sians/g' "$f"
        sed -i 's/Trros/Fo_Trros/g' "$f"
	sed -i 's/Ammex/Ba_Ammex/g' "$f"
	sed -i 's/Ammel/Ba_Ammel/g' "$f"
	sed -i 's/Asmes/Ba_Asmes/g' "$f"
	sed -i 's/Darer/Ba_Darer/g' "$f"
	sed -i 's/Gamor/Ba_Gamor/g' "$f"
	sed -i 's/Mymur/Ba_Mymur/g' "$f"
	sed -i 's/Plwal/Ba_Plwal/g' "$f"
	sed -i 's/Sigra/Ba_Sigra/g' "$f"
	sed -i 's/Thmac/Ba_Thmac/g' "$f"
	sed -i 's/Trtib/Ba_Trtib/g' "$f"
	sed -i 's/Bubuf/Ba_Bubuf/g' "$f"
	sed -i 's/Chlan/Ba_Chlan/g' "$f"
	sed -i 's/Ectel/Ba_Ectel/g' "$f"
	sed -i 's/Elmax/Ba_Elmax/g' "$f"
	sed -i 's/Ereur/Ba_Ereur/g' "$f"
	sed -i 's/Jajac/Ba_Jajac/g' "$f"
	sed -i 's/Ornil/Ba_Ornil/g' "$f"
	sed -i 's/Orcun/Ba_Orcun/g' "$f"
	sed -i 's/Ranor/Ba_Ranor/g' "$f"
	sed -i 's/Soara/Ba_Soara/g' "$f"
	sed -i 's/Sirhi/Ba_Sirhi/g' "$f"

done










