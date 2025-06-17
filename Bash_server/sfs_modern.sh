#!/bin/bash
#SBATCH -c 10
#SBATCH --mem-per-cpu 60G
#SBATCH --time=8:00:00
#SBATCH --output=/projects/mjolnir1/people/qvw641/WhoopingCrane/StairwayPlot/out/SFS_%A.log
#SBATCH --job-name SFS

module load angsd/0.940
module load winsfs/0.7.0

dir=/projects/mjolnir1/people/qvw641/WhoopingCrane/StairwayPlot/
ref=/projects/mjolnir1/people/nrv690/nrv690.hernan/GNRD_data/ref_genomes/target_species/Grus_americana_VGP/final_phased/GCF_028858705.1_bGruAme1.mat/data/GCF_028858705.1/GCF_028858705.1_bGruAme1.mat_genomic.fasta

# all modern
angsd  -b /home/qvw641/WhoopingCrane/list_bam_modernRG.txt -ref $ref -anc $ref -out ${dir}Modern_foldedSFS -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1  -C 50 -baq 0 -minMapQ 30 -minQ 20  -doCounts 1 -nThreads 10 -GL 2 -doSaf 1

winsfs ${dir}Modern_foldedSFS.saf.idx > ${dir}Modern_foldedSFS.ml
realSFS -fold 1  ${dir}Modern_foldedSFS.saf.idx -P 4 > ${dir}Modern_realfoldedSFS.ml

# wild
angsd  -b /home/qvw641/WhoopingCrane/list_bam_wildRG.txt -ref $ref -anc $ref -out ${dir}Wild_SFS -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1  -C 50 -baq 0 -minMapQ 30 -minQ 20  -doCounts 1 -nThreads 10 -GL 2 -doSaf 1

winsfs ${dir}Wild_SFS.saf.idx > ${dir}Wild_unfoldedSFS.ml
realSFS -fold 1  ${dir}Wild_SFS.saf.idx -P 4 > ${dir}Wild_foldedSFS.ml

