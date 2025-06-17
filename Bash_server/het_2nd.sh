#!/bin/bash
#SBATCH -c 10
#SBATCH --mem-per-cpu 2G
#SBATCH --time=3:00:00
#SBATCH --array=1-53%10
#SBATCH --output=/projects/erode/people/qvw641/WhoopingCrane/Het/out/Heterozygosity_2nd_%A_%a.log
#SBATCH --job-name Het_2

module load angsd/0.940
module load winsfs

dir=/projects/erode/people/qvw641/WhoopingCrane/Het/
ref=/projects/erode/people/nrv690/nrv690.hernan/GNRD_data/ref_genomes/target_species/Grus_americana_VGP/final_phased/GCF_028858705.1_bGruAme1.mat/data/GCF_028858705.1/GCF_028858705.1_bGruAme1.mat_genomic.fasta

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/WhoopingCrane/list_bam_all.txt | awk '{print $1}')
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/WhoopingCrane/Samples_all | awk '{print $1}')

zcat ${dir}${name}_rmtrans_GL2.mafs.gz  | awk '{print $1 , $2}' | sed 1d > ${dir}${name}_rmtrans.list
angsd sites index ${dir}${name}_rmtrans.list
echo "index done"

angsd  -i ${sample} -ref $ref -sites ${dir}/${name}_rmtrans.list -anc $ref -out ${dir}${name}_filtList_GL2_DP3 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1  -noTrans 1 -C 50 -baq 0 -minMapQ 30 -minQ 20  -setMinDepth 3 -doCounts 1 -nThreads 10 -GL 2 -doSaf 1

winsfs ${dir}${name}_filtList_GL2_DP3.saf.idx > ${dir}${name}_filter_gl2_dp3_win_FINAL.ml
