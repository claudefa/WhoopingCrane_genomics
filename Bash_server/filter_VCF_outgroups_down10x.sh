#!/bin/bash
#SBATCH -c 3
#SBATCH --mem-per-cpu 10G
#SBATCH --time=200:00:00
#SBATCH --output=/projects/mjolnir1/people/qvw641/WhoopingCrane/VCF_down10x/Concat/filter_outgroup_tabix.log
#SBATCH --job-name filter

# remove first fixed positions and  maxmissing 19 samples should have position covered
#vcftools --gzvcf WC_modern_down_outgroups.g.vcf.gz --maf 0.0001 --max-missing 0.4 --recode --recode-INFO-all --stdout | bgzip -c > WC_modern_down_outgroups.singletons.vcf.gz; tabix -p vcf WC_modern_down_outgroups.singletons.vcf.gz

#bcftools view -i 'F_MISSING=0' -S outgroups --exclude-types indels WC_modern_down_outgroups.singletons.vcf.gz | awk -v s=1 '{print $1"\t"$2-s"\t"$2}' > outgroups_missingness_tab.bed
#sortBed outgroups_missingness_tab.bed > outgroups_missingness_tab_sort.bed
tabix -h -R outgroups_missingness_tab_sort.bed  WC_modern_down_outgroups.singletons.vcf.gz | bgzip -c  > WC_modern_down_outgroups.singl.outmiss.vcf.gz; tabix -p vcf WC_modern_down_outgroups.singl.outmiss.vcf.gz

bcftools sort -o /projects/mjolnir1/people/qvw641/WhoopingCrane/VCF_down10x/Concat/WC_modern_down_outgroups.singl.outmiss_sort.vcf.gz -O z /projects/mjolnir1/people/qvw641/WhoopingCrane/VCF_down10x/Concat/WC_modern_down_outgroups.singl.outmiss.vcf.gz

