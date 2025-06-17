#!/bin/bash
#SBATCH -c 2
#SBATCH --mem-per-cpu 10G
#SBATCH --time=200:00:00
#SBATCH --output=/projects/mjolnir1/people/qvw641/WhoopingCrane/VCF_down10x/Concat/combine.log
#SBATCH --job-name combine

ref="/projects/mjolnir1/people/nrv690/nrv690.hernan/GNRD_data/ref_genomes/target_species/Grus_americana_VGP/final_phased/GCF_028858705.1_bGruAme1.mat/data/GCF_028858705.1/GCF_028858705.1_bGruAme1.mat_genomic.fasta"

DIR2="/projects/mjolnir1/people/qvw641/WhoopingCrane/VCF_down10x/Concat/"
DIR="/projects/mjolnir1/people/qvw641/WhoopingCrane/VCF/Concat/"

VCFs="${DIR2}/WC_modern_down_snpAD.g.vcf.gz ${DIR}/Gnigricollis_mapGruAme.g.vcf.gz ${DIR}/birdAnc314_mapGruAme.g.vcf.gz ${DIR}/birdAnc315_mapGruAme.g.vcf.gz ${DIR}/birdAnc316_mapGruAme.g.vcf.gz "

bcftools merge --force-samples --merge all $VCFs -O z -o ${DIR2}/WC_modern_down_outgroups.g.vcf.gz ;  tabix -p vcf ${DIR2}/WC_modern_down_outgroups.g.vcf.gz
