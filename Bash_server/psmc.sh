#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu 100G
#SBATCH --time=36:00:00
#SBATCH --output=/projects/mjolnir1/people/qvw641/WhoopingCrane/PSMC/test/out/psmc3.%A_%a.log
#SBATCH --job-name PSMC3

dir="/projects/mjolnir1/people/qvw641/WhoopingCrane/PSMC/test/"
/projects/mjolnir1/apps/conda/psmc-0.6.5/bin/psmc -N30 -t5 -r5 -p "4+30*2+4+6+10" -o ${dir}/WC_standard.psmc ${dir}/Crane1.diploid.psmcfa

echo "done 1st"
seq 100 | xargs -i echo /projects/mjolnir1/apps/conda/psmc-0.6.5/bin/psmc -N30 -t5 -r5 -p "4+30*2+4+6+10" -o ${dir}/boostrap/WC_{}_gf.sk.psmc ${dir}/Crane1.split.psmcfa | sh
