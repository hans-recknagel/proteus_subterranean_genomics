#!bin/bash

#Input is the path to folder with .sam bwa mem results.

FOLDER=/path/to/folder

SAM=$(ls $FOLDER/*.sam)

for s in $SAM;
do
	name=$(basename "$s" .sam)
	samtools view -bS "$s" > $FOLDER/$name.bam
	samtools sort $FOLDER/$name.bam > $FOLDER/$name.output_sorted.bam
	samtools index -c $FOLDER/$name.output_sorted.bam $FOLDER/$name.output_sorted.csi


done






