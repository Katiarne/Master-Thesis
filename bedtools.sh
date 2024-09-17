#!/bin/bash
#SBATCH --job-name=BEDTolls  # sensible name for the job
#SBATCH --mail-user=katinka.fjeld.arnesen@nmbu.no # Email me when job is done.
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END

#Collect "CDS" rows from gff
grep -w "CDS" /mnt/SCRATCH/kristenl/katinka_master/Salmo_salar.Ssal_v3.1.107.chr.gff3 > cds_only.gff3

#Load BEDTools
module load BEDTools

#Run BEDTools intersect
bedtools intersect -wa -wb -a cds_only.gff3 -b /mnt/SCRATCH/kristenl/katinka_master/Ssalv3.1_pangenome_pangenie_all_svs_aquafaang_samples.chr.vcf > intersect_otput
