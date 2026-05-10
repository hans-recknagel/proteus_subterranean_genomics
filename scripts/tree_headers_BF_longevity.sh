
files=$(ls longevity/fixedTrees/*nwk)
#cat $files

for f in $files;
do
        sed -i 's/Acbae/Fo_Acbae/g' "$f"
        sed -i 's/Cacar/Fo_Cacar/g' "$f"
        sed -i 's/Prang/Fo_Prang/g' "$f"
        sed -i 's/Elmax/Fo_Elmax/g' "$f"
        sed -i 's/Hosap/Fo_Hosap/g' "$f"
        sed -i 's/Lacha/Fo_Lacha/g' "$f"
        sed -i 's/Myluc/Fo_Myluc/g' "$f"
        sed -i 's/Hegla/Fo_Hegla/g' "$f"
        sed -i 's/Ororc/Fo_Ororc/g' "$f"
        sed -i 's/Rhtyp/Fo_Rhtyp/g' "$f"
        sed -i 's/Tecar/Fo_Tecar/g' "$f"
	sed -i 's/Ancor/Ba_Ancor/g' "$f"
	sed -i 's/Botau/Ba_Botau/g' "$f"
	sed -i 's/Capor/Ba_Capor/g' "$f"
	sed -i 's/Darer/Ba_Darer/g' "$f"
	sed -i 's/Chsab/Ba_Chsab/g' "$f"
	sed -i 's/Chpic/Ba_Chpic/g' "$f"
	sed -i 's/Cocri/Ba_Cocri/g' "$f"
	sed -i 's/Ercal/Ba_Ercal/g' "$f"
	sed -i 's/Modom/Ba_Modom/g' "$f"
	sed -i 's/Ovari/Ba_Ovari/g' "$f"
	sed -i 's/Ranor/Ba_Ranor/g' "$f"
	sed -i 's/Sasal/Ba_Sasal/g' "$f"
	sed -i 's/Scvul/Ba_Scvul/g' "$f"
	sed -i 's/Tagut/Ba_Tagut/g' "$f"
done




