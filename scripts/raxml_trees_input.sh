#!/bin/bash

#Input is the name of the analysis which the script should use for the path to the pal2nal results fot the analysis in question. It also takes an output dir. 

analysis=analysis_name

outdir=/path/to/$analysis
mkdir $outdir

dirs=$(ls -d /path/to/$analysis\_out/pal2nal_out/*mafft*)

for d in $dirs;
do
	files=$(ls "$d")

	for f in $files;
	do

		raxmlHPC-PTHREADS-SSE3 -T 10 -f a -x 123 -p 123 -N 100 -m GTRGAMMA -n "$f".nwk -w $outdir -s "$d"/"$f"

	done

done


