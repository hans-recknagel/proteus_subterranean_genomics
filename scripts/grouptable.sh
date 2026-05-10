#!/bin/bash

ANALYSIS=analysis_name

OGTABLE=/path/to/$ANALYSIS/Orthogroups.tsv

cd $ANALYSIS

mkdir orthogroups_tables_out

cd backtranslated_out

NSP=$(ls -d * | cat)

for n in $NSP;
do
	cd "$n"

	ogs=$(ls > /path/to/ogs.txt)

	name=$(cat /path/to/ogs.txt | sed 's/\.fa$//' > /path/to/name.txt)

	names=$(cat /path/to/name.txt)

	( IFS='|'; grep -E "^(${names[*]})\>" ) </path/to/$ANALYSIS/Orthogroups.tsv >/path/to/$ANALYSIS/orthogroups_tables_out/"$n"_orthogroups_table.tsv

	cd ..
done;
