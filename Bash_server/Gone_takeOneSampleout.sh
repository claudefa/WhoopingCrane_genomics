#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu 2G
#SBATCH --time=3:00:00
#SBATCH --array=1-37%10
#SBATCH --output=/projects/mjolnir1/people/qvw641/WhoopingCrane/GONe/boostrap/out/takeone_%A_%a.log
#SBATCH --job-name Bootstrap

dir="/projects/mjolnir1/people/qvw641/WhoopingCrane/GONe/"

line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${dir}Samples_modern | awk '{print $1}')

grep -v $line ${dir}/Samples_modern > ${dir}/boostrap/Sample_${line}.txt

vcftools --gzvcf /projects/mjolnir1/people/qvw641/WhoopingCrane/VCF/Concat/WC_all.vcf.gz --keep ${dir}/Samples_modern --max-missing 0.8 --maf 0.00000001 --recode --recode-INFO-all --stdout |bgzip -c > ${dir}/boostrap/WC_modern_${line}.vcf.gz ; tabix -p vcf ${dir}/boostrap/WC_modern_${line}.vcf.gz

/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --vcf ${dir}/boostrap/WC_modern_${line}.vcf.gz  --allow-extra-chr --double-id --recode --out WC_modern_${line}

while read chrom; do oldname=$(echo $chrom | awk '{print $1}'); newname=$(echo $chrom | awk '{print $2}');echo $newname;sed -i "s/$oldname/$newname/g" WC_modern_${line}.map; done < ../chromsomes_convertnames
