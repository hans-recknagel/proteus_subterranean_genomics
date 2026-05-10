#!/bin/bash

files=$(ls pigmentation_loss/fasta_out/selected_ogs/sp_only_headers/*fa)
#cat $files

for f in $files;
do

	sed -i 's/Ammex/Fo_Ammex/g' "$f"
	sed -i 's/Asmec/Fo_Asmec/g' "$f"
	sed -i 's/Hegla/Fo_Hegla/g' "$f"
	sed -i 's/Prang/Fo_Prang/g' "$f"
	sed -i 's/Luden/Fo_Luden/g' "$f"
	sed -i 's/Sirhi/Fo_Sirhi/g' "$f"
	sed -i 's/Sians/Fo_Sians/g' "$f"
	sed -i 's/Trros/Fo_Trros/g' "$f"
	sed -i 's/Ammel/Ba_Ammel/g' "$f"
	sed -i 's/Asmes/Ba_Asmes/g' "$f"
	sed -i 's/Darer/Ba_Darer/g' "$f"
	sed -i 's/Gamor/Ba_Gamor/g' "$f"
	sed -i 's/Hosap/Ba_Hosap/g' "$f"
	sed -i 's/Leocu/Ba_Leocu/g' "$f"
	sed -i 's/Mymur/Ba_Mymur/g' "$f"
	sed -i 's/Plwal/Ba_Plwal/g' "$f"
	sed -i 's/Sigra/Ba_Sigra/g' "$f"
	sed -i 's/Tagut/Ba_Tagut/g' "$f"
	sed -i 's/Thmac/Ba_Thmac/g' "$f"
	sed -i 's/Xetro/Ba_Xetro/g' "$f"
	sed -i 's/Zoviv/Ba_Zoviv/g' "$f"
	sed -i 's/Trtib/Ba_Trtib/g' "$f"
	sed -i 's/Chlan/Ba_Chlan/g' "$f"
	sed -i 's/Ornil/Ba_Ornil/g' "$f"
	sed -i 's/Scvul/Ba_Scvul/g' "$f"

done
