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

# create new directory
mkdir -p bedtools
cd bedtools

# Run BEDTools intersect with all columns and rows from eQTL file
### BRAIN ####
bedtools intersect -wa -wb -loj -a /mnt/users/katiarne/Master/vcf.bed -b /mnt/users/katiarne/Master/regions_brain.bed > intersect_output_brain.bed

bedtools intersect -wa -wb -loj -a /mnt/users/katiarne/Master/eQTLs_brain.bed -b /mnt/users/katiarne/Master/regions_brain.bed > intersect_output_eQTL_brain.bed

bedtools intersect -wa -wb -loj -a /mnt/users/katiarne/Master/eQTLs_0bp_brain.bed -b /mnt/users/katiarne/Master/gff_complete.bed > intersect_gff_brain.bed

### gonad ###
bedtools intersect -wa -wb -loj -a /mnt/users/katiarne/Master/vcf.bed -b /mnt/users/katiarne/Master/regions_gonad.bed > intersect_output_gonad.bed

bedtools intersect -wa -wb -loj -a /mnt/users/katiarne/Master/eQTLs_gonad.bed -b /mnt/users/katiarne/Master/regions_gonad.bed > intersect_output_eQTL_gonad.bed

bedtools intersect -wa -wb -loj -a /mnt/users/katiarne/Master/eQTLs_0bp_gonad.bed -b /mnt/users/katiarne/Master/gff_complete.bed > intersect_gff_gonad.bed

### liver ###
bedtools intersect -wa -wb -loj -a /mnt/users/katiarne/Master/vcf.bed -b /mnt/users/katiarne/Master/regions_liver.bed > intersect_output_liver.bed

bedtools intersect -wa -wb -loj -a /mnt/users/katiarne/Master/eQTLs_liver.bed -b /mnt/users/katiarne/Master/regions_liver.bed > intersect_output_eQTL_liver.bed

bedtools intersect -wa -wb -loj -a /mnt/users/katiarne/Master/eQTLs_0bp_liver.bed -b /mnt/users/katiarne/Master/gff_complete.bed > intersect_gff_liver.bed

### muscle ###
bedtools intersect -wa -wb -loj -a /mnt/users/katiarne/Master/vcf.bed -b /mnt/users/katiarne/Master/regions_muscle.bed > intersect_output_muscle.bed

bedtools intersect -wa -wb -loj -a /mnt/users/katiarne/Master/eQTLs_muscle.bed -b /mnt/users/katiarne/Master/regions_muscle.bed > intersect_output_eQTL_muscle.bed

bedtools intersect -wa -wb -loj -a /mnt/users/katiarne/Master/eQTLs_0bp_muscle.bed -b /mnt/users/katiarne/Master/gff_complete.bed > intersect_gff_muscle.bed
