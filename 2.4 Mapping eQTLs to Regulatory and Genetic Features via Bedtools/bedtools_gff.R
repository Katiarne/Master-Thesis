#!/bin/bash
#SBATCH --job-name=BEDTolls  # sensible name for the job
#SBATCH --mail-user=katinka.fjeld.arnesen@nmbu.no # Email me when job is done.
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mail-type=END

#Load BEDTools
module load BEDTools

# Run BEDTools intersect with all columns and rows from eQTL file
bedtools intersect -a gff.bed -b gff_gene_cordinates.bed -wa -wb  > final_gff.bed
