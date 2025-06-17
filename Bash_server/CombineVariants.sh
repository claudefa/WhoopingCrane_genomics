#!/bin/bash

#Ref genome
ref="/projects/mjolnir1/people/nrv690/nrv690.hernan/GNRD_data/ref_genomes/target_species/Grus_americana_VGP/final_phased/GCF_028858705.1_bGruAme1.mat/data/GCF_028858705.1/GCF_028858705.1_bGruAme1.mat_genomic.fasta"

#Directory
DIR="/projects/mjolnir1/people/qvw641/WhoopingCrane/VCF_down10x/"
OUTDIR=$DIR"/MergeVCFs/"
out=${OUTDIR}"out/";
qu=${OUTDIR}"qu/";

# create directories
mkdir -p $OUTDIR
mkdir -p $out
mkdir -p $qu
module unload bcftools
module load  gsl/2.5
module load bcftools/1.20

while read chrom;
do

	VCFs=""; while read line; do VCFs="$VCFs ${DIR}/$chrom/${line}_${chrom}.vcf.gz"; done  < <(cat Samples_modern_vcf_down); # Only modern downsampled
	jobName=${qu}/WC_${chrom}.combine.sh
	echo "#!/bin/bash" > $jobName
	echo "bcftools merge --force-samples --merge all $VCFs -O z -o ${OUTDIR}/WC_modern_down_${chrom}.g.vcf.gz ;  tabix -p vcf ${OUTDIR}/WC_modern_down_${chrom}.g.vcf.gz" >> $jobName	
	chmod 755 $jobName
	echo $jobName
	sbatch -c 1 --mem-per-cpu 10G --time 18:00:00 -o ${out}/Comb_${chrom}.log --job-name Comb_$chrom -- $jobName

done <   BigScaffolds_rohan.txt
