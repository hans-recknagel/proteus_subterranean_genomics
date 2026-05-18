#!/bin/bash

#Input it the genome fasta file.

seqkit sliding -g -s 1000000000 -W 1000000000 genome.fasta -o split.genome.fa.gz
