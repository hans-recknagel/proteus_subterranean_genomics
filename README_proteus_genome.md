# Proteus subterranean traits genomics

This repository accompanies:

> Recknagel et al. (in review)  
> *Phenotypic evolution of cave animals is consistently mirrored by extensive genomic convergence and remodelling*

It contains all code required to reproduce the computational analyses, statistical tests, and figure generation presented in the study.

The project investigates the genomics of subterranean adaptation and associated phenotypes across vertebrates, including:

- subterranean lifestyle adaptation
- loss of pigmentation
- eye degeneration
- longevity evolution

The study focuses primarily on the genome of the olm (*Proteus anguinus*) and comparative analyses across vertebrate taxa.

---

# Repository structure

```text
proteus_subterranean_genomics/
│
├── README.md
├── scripts/
└── data/
```

---

# Computational workflow overview

The analyses follow a modular comparative genomics pipeline:

1. Preprocessing of coding sequences
2. Orthology inference
3. Orthogroup filtering
4. CDS–protein backtranslation
5. Multiple sequence alignment
6. Codon-aware alignment reconstruction
7. Alignment cleaning and standardisation
8. Gene tree inference
9. Gene family expansion and contraction analyses
10. Selection analyses (HyPhy)
11. Genome mapping and coordinate-based analyses
12. Amino acid convergence analyses
13. Downstream statistical analyses and figure generation

---

1. Preprocessing of coding sequences

Stop codons were manually removed from CDS sequences prior to downstream analyses.

2. Orthology inference

Orthologous gene groups were inferred using OrthoFinder

Input: 
protein FASTA  files for each species
Script: 
orthofinder.sh
Output: 
orthogroups, gene trees, and orthogroup sequence sets

3. Orthogroup filtering

Orthogroups were filtered based on completeness criteria and presence across taxa.

Script: 
orthogroups_filter.R
missing_data_percentage_in_75_percent.R
Output: 
curated orthogroup sets for downstream analyses

4. CDS–protein backtranslation

Protein sequences were matched to corresponding CDS sequences to generate codon-resolved orthologs.

Scripts: 
orthogroups_filter.R
backtranslation_cleaned_files.R
backtranslation_cleaned_files_function.R
run_backtranslation_cleaned_files_function.R
run_run_backtranslation_cleaned_files_function.sh
Input:
orthogroups
protein FASTA sequences
CDS sequences (stop-codon cleaned)
Output:
codon-aware ortholog alignments

5. Multiple sequence alignment

Alignments were generated using:
MAFFT
Script:
of_result_alignment2.sh
Output:
amino acid alignments per orthogroup

6. Codon-aware alignment reconstruction

Codon alignments were generated using:
pal2nal

Script:
pal2nal_back_script2.sh
Output:
codon-resolved alignments for downstream analyses

7. Alignment cleaning
FASTA headers were standardised and cleaned.

Script:
clean_double_names_fa.sh

8. Gene tree inference (per orthogroup)

Maximum likelihood gene trees were inferred using:
RAxML

Script:
raxml_trees_input.sh
Output:
best-scoring ML trees per orthogroup

These trees are reused across selection, convergence, and comparative analyses.

9. Gene family expansion and contraction analyses

Gene family size evolution was inferred from orthogroup counts.

Scripts:
contracted_expanded_gene_groups_filtering.R
Output:
expanded and contracted gene families
summary statistics of gene family evolution

10. Selection analyses 

Foreground and background lineages were defined for hypothesis-driven analyses:

Scripts:
foreground_cave_species.sh
foreground_eye_loss.sh
foreground_longevity.sh
foreground_pigmentation_loss.sh
mark_*_relax.sh

Selection analyses were performed using:
HyPhy
aBSREL (episodic selection)
RELAX (selection intensity shifts)

Script:
hyphy_*_results.R
hyphy_absrel.sh
hyphy_relax.sh
RELAX
Post-processing:
rerun_empty_relax_results.sh
rerun_nan_relax_results.sh

11. Genome mapping and coordinate-based analyses

Read mapping and coordinate-based analyses were performed using:
BWA
SAMtools

Scripts:
filter_for_prang.sh
bwa_mem.sh
sam_to_bam.sh
bwa_mem_og_to_genome_alignments_to_table.R
split_by_GB.sh
gff_for_split_genome.R

12. Amino acid convergence analyses

Pipeline:
FASTA header standardisation
fasta_headers_all_OGs.sh
Foreground/background annotation
fasta_headers_BF_*
Format conversion
seqret.sh
Tree preparation
copy_trees.sh
tree_headers_BF_*
Convergence inference
TDG09 pipeline
tdg09_timeout.sh
Result parsing
tdg09_results.R

13. Downstream statistical analyses and figure generation

All statistical analyses and figure/table generation were performed in R.

Script:
Proteus_genome_analysis_cleaned_R.R
Input:
processed outputs from all upstream analyses
Output:
all figures and statistical results presented in the manuscript stored in the results/ directory

---

Reproducibility

All analyses were conducted using relative file paths anchored to the repository root to ensure portability across systems.

The workflow was executed on Linux-based HPC systems. Computationally intensive steps (e.g. HyPhy analyses and convergence inference) were parallelised where appropriate.
Software requirements

---

Software requirements

All required R libraries are specified at the beginning of each R script.

---

Contact

Hans Recknagel
University of Trier
recknagel@uni-trier.de

Luka Močivnik
Biotechnical Faculty, University of Ljubljana
luka.mocivnik@bf.uni-lj.si

