#!/bin/bash

#Input are the properly prepared AA fasta files. 

orthofinder -t 40 \
	-a 8 \
	-M msa \
	-S blast \
	-A mafft \
	-n OrthoFinder_cave_species \
	-f /abs/path/to/pep/files;
