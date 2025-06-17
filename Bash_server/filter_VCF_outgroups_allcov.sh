#!/bin/bash
#SBATCH -c 3
#SBATCH --mem-per-cpu 10G
#SBATCH --time=30:00:00
#SBATCH --output=/projects/mjolnir1/people/qvw641/WhoopingCrane/VCF/Concat/filter_outgroup_tabix.log
#SBATCH --job-name filter

#bcftools view -m2 -M2 -v snps WC_all.g.vcf.gz | bgzip -c > WC_all_filter.g.vcf.gz
#/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --vcf WC_all_filter.g.vcf.gz  --allow-extra-chr --double-id --maf 0.05 --geno 0.2 --recode --out WC_all_filter_gvcf

#/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --file WC_all_filter_gvcf  --pca --out WC_all_filter_gvcf --allow-extra-chr --double-id

########3

# remove first fixed positions and  maxmissing 19 samples should have position covered
#vcftools --gzvcf WC_outgroups.g.vcf.gz --maf 0.0001 --max-missing 0.4 --recode --recode-INFO-all --stdout | bgzip -c > WC_outgroups.singletons.vcf.gz; tabix -p vcf WC_outgroups.singletons.vcf.gz

#tabix -h -R outgroups_missingness_tab.bed  WC_outgroups.singletons.vcf.gz | bgzip -c  > WC_outgroups.singl.outmiss2.vcf.gz; tabix -p vcf WC_outgroups.singl.outmiss2.vcf.gz


vcftools --gzvcf WC_outgroups.g.vcf.gz --keep Samples_modern_outgroups --maf 0.0001 --max-missing 0.4 --recode --recode-INFO-all --stdout | bgzip -c > WC_modern_outgroups.singletons.vcf.gz; tabix -p vcf WC_modern_outgroups.singletons.vcf.gz

