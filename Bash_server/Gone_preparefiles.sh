vcftools --gzvcf /projects/mjolnir1/people/qvw641/WhoopingCrane/VCF/Concat/WC_all.vcf.gz --keep Samples_wild --max-missing 0.8 --maf 0.0000001 --recode --recode-INFO-all --stdout |bgzip -c > WC_wild.vcf.gz ; tabix -p vcf WC_wild.vcf.gz

vcftools --gzvcf /projects/mjolnir1/people/qvw641/WhoopingCrane/VCF/Concat/WC_all.vcf.gz --keep Samples_captive_earlyfounders --max-missing 0.8 --maf 0.0000001 --recode --recode-INFO-all --stdout |bgzip -c > WC_captive_earlyfounders.vcf.gz ; tabix -p vcf WC_captive_earlyfounders.vcf.gz

vcftools --gzvcf /projects/mjolnir1/people/qvw641/WhoopingCrane/VCF/Concat/WC_all.vcf.gz --keep Samples_captive --max-missing 0.8 --maf 0.00000001 --recode --recode-INFO-all --stdout |bgzip -c > WC_captive.vcf.gz ; tabix -p vcf WC_captive.vcf.gz


vcftools --gzvcf /projects/mjolnir1/people/qvw641/WhoopingCrane/VCF/Concat/WC_all.vcf.gz --keep Samples_modern --max-missing 0.8 --maf 0.00000001 --recode --recode-INFO-all --stdout |bgzip -c > WC_modern.vcf.gz ; tabix -p vcf WC_modern.vcf.gz


/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --vcf WC_wild.vcf.gz  --allow-extra-chr --double-id --recode --out WC_wild
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --vcf WC_captive_earlyfounders.vcf.gz  --allow-extra-chr --double-id --recode --out WC_captive_earlyfounders
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --vcf WC_captive.vcf.gz  --allow-extra-chr --double-id --recode --out WC_captive
/projects/mjolnir1/apps/conda/plink-1.90b6.21/bin/plink --vcf WC_modern.vcf.gz  --allow-extra-chr --double-id --recode --out WC_modern


# change chromosome names 

while read line; do oldname=$(echo $line | awk '{print $1}'); newname=$(echo $line | awk '{print $2}');echo $newname;sed -i "s/$oldname/$newname/g" WC_wild.map; done < chromsomes_convertnames
while read line; do oldname=$(echo $line | awk '{print $1}'); newname=$(echo $line | awk '{print $2}');echo $newname;sed -i "s/$oldname/$newname/g" WC_captive.map; done < chromsomes_convertnames
while read line; do oldname=$(echo $line | awk '{print $1}'); newname=$(echo $line | awk '{print $2}');echo $newname;sed -i "s/$oldname/$newname/g" WC_captive_earlyfounders.map; done < chromsomes_convertnames
