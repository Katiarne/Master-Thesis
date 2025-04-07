#!/bin/bash
#SBATCH --job-name=BEDTools  # sensible name for the job
#SBATCH --mail-user=katinka.fjeld.arnesen@nmbu.no # Email me when job is done.
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mail-type=END

#Load BEDTools
module load BEDTools

# Run BEDTools intersect with all columns and rows from eQTL file
bedtools intersect -wa -wb -loj -a vcf.bed -b gff_complete.bed > intersect_output_gff.bed

bedtools intersect -wa -wb -loj -a vcf.bed -b regions.bed > intersect_output_all.bed

bedtools intersect -wa -wb -loj -a eQTLs.bed -b regions.bed > intersect_output.bed

bedtools intersect -wa -wb -loj -a eQTLs_0bp.bed -b gff_complete.bed > intersect_gff_complete.bed
