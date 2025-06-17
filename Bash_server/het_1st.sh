#!/bin/bash
#SBATCH -c 2
#SBATCH --mem-per-cpu 10G
#SBATCH --time=8:00:00
#SBATCH --array=1-53%20
#SBATCH --output=/projects/erode/people/qvw641/WhoopingCrane/Het/out/Heterozygosity_1st_%A_%a.log
#SBATCH --job-name Het_1st

module load angsd/0.940
dir=/projects/erode/people/qvw641/WhoopingCrane/Het/
ref=/projects/erode/people/nrv690/nrv690.hernan/GNRD_data/ref_genomes/target_species/Grus_americana_VGP/final_phased/GCF_028858705.1_bGruAme1.mat/data/GCF_028858705.1/GCF_028858705.1_bGruAme1.mat_genomic.fasta

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/WhoopingCrane/list_bam_all.txt | awk '{print $1}')
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/WhoopingCrane/Samples_all | awk '{print $1}')
cov=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/WhoopingCrane/Samples_all | awk '{print $2}')

angsd  -i $sample -ref $ref -rf /home/qvw641/WhoopingCrane/BigScaffolds.txt -anc $ref  -out ${dir}${name}_rmtrans_GL2 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -rmTrans 1 -trim 0 -C 50 -baq 0 -minMapQ 20 -minQ 20 -setMinDepth 3 -setMaxDepth $cov -doCounts 1 -doMajorMinor 1 -GL 2 -doGlf 2 -doMaf 2
