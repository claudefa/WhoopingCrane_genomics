#!/bin/bash
#SBATCH -c 16
#SBATCH --mem-per-cpu 2G
#SBATCH --time=48:00:00
#SBATCH --array=1-33%10
#SBATCH --output=/projects/mjolnir1/people/qvw641/WhoopingCrane/ROH_down10x/out/Rohan_%A_%a.log
#SBATCH --job-name Rohan

source /home/qvw641/bin/rohan_env/bin/activate
dir=/projects/mjolnir1/people/qvw641/WhoopingCrane/ROH_down10x/

ref="/projects/mjolnir1/people/nrv690/nrv690.hernan/GNRD_data/ref_genomes/target_species/Grus_americana_VGP/final_phased/GCF_028858705.1_bGruAme1.mat/data/GCF_028858705.1/GCF_028858705.1_bGruAme1.mat_genomic.fasta"

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/WhoopingCrane/Down_bams_10x | awk '{print $1}')
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/WhoopingCrane/DownSample_10x | awk '{print $1}')

# modern without damage
/home/qvw641/bin/rohan/src/rohan -t 16  --size 1000000 --rohmu 2e-5 --auto /home/qvw641/WhoopingCrane/BigScaffolds_rohan.txt -o ${dir}/${name}_modern_2e5_1Mb_nochrZ_down10x $ref ${sample}


