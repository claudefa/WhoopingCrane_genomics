#!/bin/bash
#SBATCH -c 12
#SBATCH --mem-per-cpu 10G
#SBATCH --time=12:00:00  
##SBATCH --array=1-37%1
#SBATCH --output=/projects/mjolnir1/people/qvw641/WhoopingCrane/VCF_down10x/out/Repeat_snpad1.%A_%a.log
#SBATCH --job-name snpad


# Set directories
DIR="/projects/mjolnir1/people/qvw641/WhoopingCrane/"
OUTDIR=${DIR}"VCF_down10x/"

sample=$2
bam=$(grep -w $sample /home/qvw641/WhoopingCrane/Samples_down10x_calling | awk '{print $2}')
chrom=$1




vcf=${OUTDIR}/${chrom}/${sample}_${chrom}.vcf
#REF
ASSEMBLY="/projects/mjolnir1/people/nrv690/nrv690.hernan/GNRD_data/ref_genomes/target_species/Grus_americana_VGP/final_phased/GCF_028858705.1_bGruAme1.mat/data/GCF_028858705.1/GCF_028858705.1_bGruAme1.mat_genomic.fasta"
module load jdk/1.8.0_291
module load picard/2.27.5
module load snpAD/0.3.10

snpad=${OUTDIR}/${chrom}/${sample}_${chrom}.snpAD
mkdir -p ${OUTDIR}/${chrom}/

Bam2snpAD -r $chrom -f $ASSEMBLY -Q 30 /projects/mjolnir1/people/qvw641/WhoopingCrane/BAMs_10x/$bam  > $snpad
snpAD -c 12 -o ${OUTDIR}/${chrom}/${sample}.priors.txt -O ${OUTDIR}/${chrom}/${sample}.errors.txt $snpad
echo "snpAD 1st done"
snpADCall -N 1 -e ${OUTDIR}/${chrom}/${sample}.errors.txt -p "`cat ${OUTDIR}/${chrom}/${sample}.priors.txt`" $snpad > $vcf
rm ${vcf}.gz
bgzip $vcf; tabix -p vcf ${vcf}.gz
