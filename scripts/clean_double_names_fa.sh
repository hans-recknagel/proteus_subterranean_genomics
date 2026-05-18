#!/bin/bash

#Input are pal2nal results for a particular analysis.

ANALYSIS=analysis_name

FOLDERS=$(ls -d $ANALYSIS\_out/pal2nal_out/*)
for f in $FOLDERS;
do
	INFILES=$(ls "$f"/*)
		for i in $INFILES;
		do
			sed -i 's/^>Nefor_/>/' "$i"
			sed -i 's/^>Menas_/>/' "$i"
			sed -i 's/^>Pipip_/>/' "$i"
		done

done
