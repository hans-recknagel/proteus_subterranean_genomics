#!/bin/bash

analysis=analysis_name

#OGS=$analysis\_ogs_all_OGonly.txt
OGS=$analysis\_failed_jobs.log

treedir=$analysis/fixedTrees
phyldir=$analysis/phylip
outdir=$analysis/tdg09_out_timeout
#LOG_FILE=$analysis\_done_jobs_rerun.log
LOG_FILE=$analysis\_done_jobs.log
#FAILED_LOG=$analysis\_failed_jobs.log
FAILED_LOG=$analysis\_failed_jobs_rerun.log

rm $LOG_FILE
rm $FAILED_LOG

touch "$LOG_FILE"
touch "$FAILED_LOG"

mkdir $outdir

og=$(cat $OGS)

for o in $og;
do
	
	if grep -Fxq "$o" "$LOG_FILE"; then
        echo "Skipping $filename (already completed)"
        continue
    fi
	
	if timeout 300 java -cp ~/software/tdg09/tdg09-1.1.2/dist/tdg09.jar tdg09.Analyse -alignment $phyldir/aligned_mafft_"$o"_filtered_sp_only.phylip -groups Fo Ba -threads 4 -tree $treedir/RAxML_bestTree."$o".fa_aligned_mafft_pal2nal_out.fa.nwk > $outdir/"$o".tdg09.out;
  then
        echo "$o" >> "$LOG_FILE"
        echo "Completed $o"
    else
        echo "Failed or timed out: $filename"
	echo "$o" >> "$FAILED_LOG"

    fi
	
done;

echo "All files processed or attempted."

