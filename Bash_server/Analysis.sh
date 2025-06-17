module load angsd/0.940
ref="/projects/mjolnir1/people/nrv690/nrv690.hernan/GNRD_data/ref_genomes/target_species/Grus_americana_VGP/final_phased/GCF_028858705.1_bGruAme1.mat/data/GCF_028858705.1/GCF_028858705.1_bGruAme1.mat_genomic.fasta"


#Relatedness only historical and only modern ---------
#Historical
dir="/projects/mjolnir1/people/qvw641/WhoopingCrane/Inbreeding/"
jobName=$dir/out/angsd_inbreeding_tvonly_onlyHistorical.sh
echo '#!/bin/bash' > $jobName
echo "angsd -bam /home/qvw641/WhoopingCrane/list_bam_historical_newseq.txt -ref $ref -uniqueOnly 1 -remove_bads 1 -noTrans 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 0 -skipTriallelic 1 -gl 2 -minMapQ 30 -nThreads 10 -doGlf 3 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -P 4 -out ${dir}/WC_all_noTrans_historical" >> $jobName
sbatch -c 4 --mem-per-cpu 10G --time 48:00:00 -o ${dir}/out/ANGSD_WC_rel_notrans_historical.log --job-name Rel_noTrans -- $jobName

zcat WC_all_noTrans_historical.mafs.gz | cut -f6 | sed 1d > WC_historical_noTrans.freq
jobName=$dir/out/NGSrelate_rel_historical.sh
echo '#!/bin/bash' > $jobName
echo "/projects/mjolnir1/apps/ngsRelate/ngsRelate  -g ${dir}WC_all_noTrans_historical.glf.gz -n 16 -f ${dir}WC_historical_noTrans.freq -O ${dir}WC_historical_noTrans_rel.ml" >> $jobName
sbatch -c 1 --mem-per-cpu 2G --time 12:00:00 -o ${dir}/out/NGSrelate_WC_rel.log --job-name Rel_WC -- $jobName


#Modern
jobName=$dir/out/angsd_inbreeding_allpos_onlyModern.sh
echo '#!/bin/bash' > $jobName
echo "angsd -bam /home/qvw641/WhoopingCrane/list_bam_modern.txt -ref $ref -uniqueOnly 1 -remove_bads 1 -noTrans 0 -only_proper_pairs 1 -trim 0 -C 50 -baq 0 -skipTriallelic 1 -gl 2 -minMapQ 30 -nThreads 10 -doGlf 3 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -P 4 -out ${dir}/WC_all_allPos_modern" >> $jobName
sbatch -c 4 --mem-per-cpu 10G --time 48:00:00 -o ${dir}/out/ANGSD_WC_rel_all_modern.log --job-name Rel_all -- $jobName

zcat WC_all_allPos_modern.mafs.gz | cut -f6 | sed 1d > WC_modern_allPos.freq
jobName=$dir/out/NGSrelate_rel_modern.sh
echo '#!/bin/bash' > $jobName
echo "/projects/mjolnir1/apps/ngsRelate/ngsRelate  -g ${dir}WC_all_allPos_modern.glf.gz -n 37 -f ${dir}WC_modern_allPos.freq -O ${dir}WC_modern_allPos_rel.ml" >> $jobName
sbatch -c 1 --mem-per-cpu 2G --time 12:00:00 -o ${dir}/out/NGSrelate_WC_rel.log --job-name RelC_modern -- $jobName


# Heterozygosity ---------------------
# Calculate heterozygosity from SFS
# first pre-run ANGSD with rmTrans 1 to extract the positions

module load angsd/0.940
het_1st.sh
het_2nd.sh

# VCF snpAD ---------------
bash run_vcf.sh # runs snpAD with historical (snpAD_historical.sh) and modern (snpAD_modern.sh)
bash CombineVariants.sh # combine variants per chromosome
bash filtervcf.sh # basic filters
bash concatVCFs_filter.sh # or for gvcfs concatgVCFs_filter.sh 

