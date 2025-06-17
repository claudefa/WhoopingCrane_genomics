#!/bin/bash
#SBATCH -c 4
#SBATCH --mem-per-cpu 10G
#SBATCH --time=24:00:00
#SBATCH --array=1-37
#SBATCH --output=/projects/mjolnir1/people/qvw641/WhoopingCrane/GONe/AllModern_3.42_new/boostrap/out_job/gone_bs_%A_%a.log
#SBATCH --job-name gone_bs

dir="/projects/mjolnir1/people/qvw641/WhoopingCrane/GONe/AllModern_3.42_new/"

line=$(sed -n "$SLURM_ARRAY_TASK_ID"p /projects/mjolnir1/people/qvw641/WhoopingCrane/GONe/Samples_modern | awk '{print $1}')

mkdir $line
cd ${dir}/boostrap/${line}
cp ${dir}/boostrap/WC_modern_${line}.map ${dir}/boostrap/${line}
cp ${dir}/boostrap/WC_modern_${line}.ped ${dir}/boostrap/${line}

bash ${dir}/script_GONE.sh WC_modern_$line
