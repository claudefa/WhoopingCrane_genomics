#!/bin/bash

DIR="/projects/mjolnir1/people/qvw641/WhoopingCrane/"
IN="/projects/mjolnir1/people/nrv690/nrv690.hernan/GNRD_analysis/cranes/historical/mapping_manual/USCS_method_final_ref_GCF_028858705.1_bGruAme1.mat/0_merged_libraries/"
IN2="/projects/mjolnir1/people/nrv690/nrv690.hernan/GNRD_analysis/cranes/modern/final_ref_GCF_028858705.1_bGruAme1.mat/"
OUTDIR=${DIR}"/Coverage/"

qu=$OUTDIR"qu"
out=${OUTDIR}"out"

mkdir -p $OUTDIR
mkdir -p $qu
mkdir -p $out
module load mosdepth/0.3.3

# send array
array=$qu/Coverage.sh
echo "#!/bin/bash" > $array
echo "#" >> $array
echo "#SBATCH -c 1" >> $array
echo "#SBATCH --mem-per-cpu 3G" >> $array
echo "#SBATCH --time=12:00:00"	>> $array
echo "#SBATCH --array=12%10" >> $array
echo "#SBATCH --output=${out}/Mosdepth_historical.%A_%a.log" >> $array
echo "#SBATCH --job-name MD_CF" >> $array

echo "module load mosdepth/0.3.3" >> $array
echo "id=\$(sed -n \"\$SLURM_ARRAY_TASK_ID\"p /home/qvw641/WhoopingCrane/Samples_historical | awk '{print \$1}')" >> $array
#echo "id=\$(sed -n \"\$SLURM_ARRAY_TASK_ID\"p /home/qvw641/WhoopingCrane/Samples_modern | awk '{print \$1}')" >> $array
echo "echo \$id" >> $array

echo "mosdepth -n --by 100000 ${OUTDIR}/\$id ${IN}/\${id}_merged.bam" >> $array
#echo "mosdepth -n --by 100000 ${OUTDIR}/\$id ${IN2}/\${id}.craneVGP_GCF_028858705.1_bGruAme1.mat.bam" >> $array

sbatch $array
