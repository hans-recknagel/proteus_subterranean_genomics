#!/bin/bash

INFOLDER=/path/to/filtered/and/backtranslated/data
#OUTFOLDER=data
#aligner is MUSCLE, MAFFT or BOTH
ALIGNER=MAFFT

cd $INFOLDER

DIRS=$(ls | cat)
#echo $DIRS

for d in $DIRS;
do
	cd "$d"
	INDIRS=$(ls | cat)
#	echo $INDIRS
	for i in $INDIRS;
	do
#		if grep -qw "MUSCLE" $ALIGNER
		if [[ $ALIGNER == "MUSCLE" ]];
		then
			mkdir "$i"_aligned_muscle
			mkdir "$i"_aligned_muscle_stable
#		elif grep -qw "MAFFT" $ALIGNER
		elif [[ $ALIGNER == "MAFFT" ]];
		then
			mkdir "$i"_aligned_mafft
		elif [[ $ALIGNER == "BOTH"  ]]
		then
			mkdir "$i"_aligned_muscle
			mkdir "$i"_aligned_muscle_stable
			mkdir "$i"_aligned_mafft
		else
			echo "Wrong aligner input."
			exit 1
		fi;
		cd "$i"
		FILES=$(ls | cat)
#		echo $FILES
		for f in $FILES;
		do
#			echo "$f"
			if [[ $ALIGNER == 'MUSCLE' ]]
			then
				/data/SOFTWARE/MUSCLE/muscle5.1 -align "$f" -output ../"$i"_aligned_muscle/aligned_muscle_"$f"
				python /data/home/lukam/software/stable.py "$f" ../"$i"_aligned_muscle/aligned_muscle_"$f" > ../"$i"_aligned_muscle_stable/aligned_muscle_stable_"$f"
			elif [[ $ALIGNER == 'MAFFT' ]]
			then
				mafft --thread 20 "$f" > ../"$i"_aligned_mafft/aligned_mafft_"$f"
			elif [[ $ALIGNER == "BOTH"  ]]
			then
				/data/SOFTWARE/MUSCLE/muscle5.1 -align "$f" -output ../"$i"_aligned_muscle/aligned_muscle_"$f"
				python /data/home/lukam/software/stable.py "$f" ../"$i"_aligned_muscle/aligned_muscle_"$f" > ../"$i"_aligned_muscle_stable/aligned_muscle_stable_"$f"
				mafft --thread 20 "$f" > ../"$i"_aligned_mafft/aligned_mafft_"$f"

			else
				echo "Wrong aligner input."
			fi;
		done;
	cd ..

	done;
	cd ..
done;
