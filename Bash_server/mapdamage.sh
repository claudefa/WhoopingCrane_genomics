#!/bin/bash

DIR="/projects/mjolnir1/people/qvw641/WhoopingCrane/"
IN="/projects/mjolnir1/people/nrv690/nrv690.hernan/GNRD_analysis/cranes/historical/mapping_manual/USCS_method_final_ref_GCF_028858705.1_bGruAme1.mat/0_merged_libraries/"
IN2="/projects/mjolnir1/people/nrv690/nrv690.hernan/GNRD_analysis/cranes/modern/final_ref_GCF_028858705.1_bGruAme1.mat/"
OUTDIR=${DIR}"/MapDamage/"

qu=$OUTDIR"qu"
out=${OUTDIR}"out"

mkdir -p $OUTDIR
mkdir -p $qu
mkdir -p $out
module load mapdamage2/

# send array
array=$qu/Mapdamage.sh
echo "#!/bin/bash" > $array
echo "#" >> $array
echo "#SBATCH -c 1" >> $array
echo "#SBATCH --mem-per-cpu 10G" >> $array
echo "#SBATCH --time=12:00:00"	>> $array
echo "#SBATCH --array=12%10" >> $array
echo "#SBATCH --output=${out}/Damage_modern.%A_%a.log" >> $array
echo "#SBATCH --job-name MapDamage" >> $array

echo "module load mapdamage2/" >> $array
echo "id=\$(sed -n \"\$SLURM_ARRAY_TASK_ID\"p /home/qvw641/WhoopingCrane/Samples_historical | awk '{print \$1}')" >> $array
#echo "id=\$(sed -n \"\$SLURM_ARRAY_TASK_ID\"p /home/qvw641/WhoopingCrane/Samples_modern | awk '{print \$1}')" >> $array
echo "ref=/projects/mjolnir1/people/nrv690/nrv690.hernan/GNRD_data/ref_genomes/target_species/Grus_americana_VGP/final_phased/GCF_028858705.1_bGruAme1.mat/data/GCF_028858705.1/GCF_028858705.1_bGruAme1.mat_genomic.fasta" >> $array
echo "mapDamage -i ${IN}/\${id}_merged.bam -r \$ref" >> $array
#echo "mapDamage -i ${IN2}/\${id}.craneVGP_GCF_028858705.1_bGruAme1.mat.bam -r \$ref" >> $array

sbatch $array
