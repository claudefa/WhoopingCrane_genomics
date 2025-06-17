#!/bin/bash

#Script to filter variants per chromosomes 

#Directory
DIR="/projects/mjolnir1/people/qvw641/WhoopingCrane/VCF/"
DIR="/projects/mjolnir1/people/qvw641/WhoopingCrane/VCF_down10x/"
INDIR=$DIR"/MergeVCFs/"
OUTDIR=$DIR"/FilteredVCF/"
out=${OUTDIR}"out/";
qu=${OUTDIR}"qu/";

# create directories
mkdir -p $OUTDIR
mkdir -p $out
mkdir -p $qu


module load bcftools
while read line
do
	echo $line
	jobName=$qu/Filter_${line}.sh
        
	echo '#!/bin/bash' > $jobName
	#echo "module load bcftools; bcftools filter -e 'FORMAT/DP<4 && FORMAT/GQ<30 && FORMAT/DP>50' ${INDIR}/WC_all_${line}.g.vcf.gz | bgzip -c > ${OUTDIR}/WC_all_${line}_filter_dp4.g.vcf.gz; tabix -p vcf ${OUTDIR}/WC_all_${line}_filter_dp4.g.vcf.gz" >> $jobName 
        echo "module load bcftools; bcftools filter -e 'FORMAT/DP<4 && FORMAT/GQ<30 && FORMAT/DP>50' ${INDIR}/WC_modern_down_${line}.g.vcf.gz | bgzip -c > ${OUTDIR}/WC_all_${line}_filter.g.vcf.gz; tabix -p vcf ${OUTDIR}/WC_modern_down_${line}_filter.g.vcf.gz" >> $jobName 
	chmod 755 $jobName
	sbatch -c 1 --mem-per-cpu 10G --time 12:00:00 -o ${out}/Filt_${line}_all.log --job-name filt$line -- $jobName	
done < BigScaffolds_rohan.txt
