#!bin/bash

FOLDER=bwa_mem_out

SAM=$(ls $FOLDER/*.sam)

for s in $SAM;
do
	name=$(basename "$s" .sam)
	samtools view -bS "$s" > $FOLDER/$name.bam
	samtools sort $FOLDER/$name.bam > $FOLDER/$name.output_sorted.bam
	samtools index -c $FOLDER/$name.output_sorted.bam $FOLDER/$name.output_sorted.csi


done






