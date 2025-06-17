#!/bin/bash
#SBATCH -c 2
#SBATCH --mem-per-cpu 10G
#SBATCH --time=200:00:00
#SBATCH --output=/projects/mjolnir1/people/qvw641/WhoopingCrane/VCF/Concat/combine.log
#SBATCH --job-name combine

bcftools view -m2 -M2 -v snps WC_all.g.vcf.gz | bgzip -c > WC_all_filter.g.vcf.gz


ref="/projects/mjolnir1/people/nrv690/nrv690.hernan/GNRD_data/ref_genomes/target_species/Grus_americana_VGP/final_phased/GCF_028858705.1_bGruAme1.mat/data/GCF_028858705.1/GCF_028858705.1_bGruAme1.mat_genomic.fasta"

DIR="/projects/mjolnir1/people/qvw641/WhoopingCrane/VCF/Concat/"

module load bcftools

VCFs="${DIR}/WC_all.g.vcf.gz ${DIR}/Gnigricollis_mapGruAme.g.vcf.gz ${DIR}/birdAnc314_mapGruAme.g.vcf.gz ${DIR}/birdAnc315_mapGruAme.g.vcf.gz ${DIR}/birdAnc316_mapGruAme.g.vcf.gz "

bcftools merge --force-samples --merge all $VCFs -O z -o ${DIR}/WC_outgroups.g.vcf.gz ;  tabix -p vcf ${DIR}/WC_outgroups.g.vcf.gz
