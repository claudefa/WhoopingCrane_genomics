#!/bin/bash
#SBATCH -c 12
#SBATCH --mem-per-cpu 10G
#SBATCH --time=34:00:00  
#SBATCH --array=1-7%2
#SBATCH --output=/projects/mjolnir1/people/qvw641/WhoopingCrane/VCF/out/snpad_historical.%A_%a.log
#SBATCH --job-name snpad


# Set directories
DIR="/projects/mjolnir1/people/qvw641/WhoopingCrane/"
OUTDIR=${DIR}"VCF/"

bam=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/WhoopingCrane/list_bam_historical_newseq_highCov.txt | awk '{print $1}')
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/WhoopingCrane/Samples_historical_highCov | awk '{print $1}')
chrom=$1
vcf=${OUTDIR}/${chrom}/${sample}_${chrom}.vcf
#REF
ASSEMBLY="/projects/mjolnir1/people/nrv690/nrv690.hernan/GNRD_data/ref_genomes/target_species/Grus_americana_VGP/final_phased/GCF_028858705.1_bGruAme1.mat/data/GCF_028858705.1/GCF_028858705.1_bGruAme1.mat_genomic.fasta"
module load jdk/1.8.0_291
module load picard/2.27.5
module load snpAD/0.3.10

snpad=${OUTDIR}/${chrom}/${sample}_${chrom}.snpAD
mkdir -p ${OUTDIR}/${chrom}/

Bam2snpAD -r $chrom -f $ASSEMBLY -Q 30 $bam  > $snpad
snpAD -c 12 -o ${OUTDIR}/${chrom}/${sample}.priors.txt -O ${OUTDIR}/${chrom}/${sample}.errors.txt $snpad
echo "snpAD 1st done"
snpADCall -N 1 -e ${OUTDIR}/${chrom}/${sample}.errors.txt -p "`cat ${OUTDIR}/${chrom}/${sample}.priors.txt`" $snpad > $vcf
bgzip $vcf; tabix -p vcf ${vcf}.gz
