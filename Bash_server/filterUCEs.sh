#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu 2G
#SBATCH --time=48:00:00
#SBATCH --output=/projects/mjolnir1/people/qvw641/WhoopingCrane/Annotation/out/Filter_intersect_%A.log
#SBATCH --job-name FilterUCEs

ref="/projects/mjolnir1/people/nrv690/nrv690.hernan/GNRD_data/ref_genomes/target_species/Grus_americana_VGP/final_phased/GCF_028858705.1_bGruAme1.mat/data/GCF_028858705.1/GCF_028858705.1_bGruAme1.mat_genomic.fasta"

module load openjdk/17.0.8
 module load python
 module load gatk
gatk SelectVariants -R $ref -L /projects/mjolnir1/people/qvw641/WhoopingCrane/Annotation/UCE_GruAme_regions.bed -V  /projects/mjolnir1/people/qvw641/WhoopingCrane/Annotation/WC_outgroups.singl.outmiss.vcf.gz  -O /projects/mjolnir1/people/qvw641/WhoopingCrane/Annotation/WC_outgroup_UCE.ann.vcf.gz; 
tabix -p vcf /projects/mjolnir1/people/qvw641/WhoopingCrane/Annotation/WC_outgroup_UCE.ann.vcf.gz

