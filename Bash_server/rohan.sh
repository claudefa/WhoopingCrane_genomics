#!/bin/bash
#SBATCH -c 16
#SBATCH --mem-per-cpu 2G
#SBATCH --time=8:00:00
#SBATCH --array=1-53%10
#SBATCH --output=/projects/mjolnir1/people/qvw641/WhoopingCrane/ROH/out/Rohan_%A_%a.log
#SBATCH --job-name Rohan

source /home/qvw641/bin/rohan_env/bin/activate
dir=/projects/mjolnir1/people/qvw641/WhoopingCrane/ROH/

ref="/projects/mjolnir1/people/nrv690/nrv690.hernan/GNRD_data/ref_genomes/target_species/Grus_americana_VGP/final_phased/GCF_028858705.1_bGruAme1.mat/data/GCF_028858705.1/GCF_028858705.1_bGruAme1.mat_genomic.fasta"

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/WhoopingCrane/list_bam_all.txt | awk '{print $1}')
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/WhoopingCrane/Samples_all | awk '{print $1}')

# Historical
/home/qvw641/bin/rohan/src/estimateDamage.pl --length 50 --threads 16 -o $dir/${name} $ref ${sample}
/home/qvw641/bin/rohan/src/rohan -t 16  --size 500000 --rohmu 2e-5 --deam5p ${dir}/${name}.5p.prof --deam3p ${dir}/${name}.3p.prof --auto /home/qvw641/WhoopingCrane/BigScaffolds_rohan.txt -o ${dir}/${name}_aDNA_2e5_500Kb_nochrZ $ref ${sample}


# modern without damage
/home/qvw641/bin/rohan/src/rohan -t 16  --size 1000000 --rohmu 2e-5 --auto /home/qvw641/WhoopingCrane/BigScaffolds_rohan.txt -o ${dir}/${name}_modern_2e5_1Mb_nochrZ $ref ${sample}


