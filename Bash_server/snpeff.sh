# Prepare the genome, copy fasta file, gff and cds to /home/qvw641/bin/snpEff/data/GrusAme and change the name

mv genomic.gff genes.gff
mv cds_from_genomic.fna sequences.fa
mv GCF_028858705.1_bGruAme1.mat_genomic.fasta GrusAme.fa


# run 
java -jar snpEff.jar  build -gff3 -v GrusAme 

./annotatevcf.sh

# Extract different profiles
zgrep -E "#|HIGH" WC_outgroups.singl.outmiss.vcf.gz | bcftools query -f '%CHROM %POS  %REF  %ALT [\t%GT:%DP]\n' > high.txt
zgrep -E "#|MODERATE" WC_outgroups.singl.outmiss.vcf.gz | bcftools query -f '%CHROM %POS  %REF  %ALT [\t%GT:%DP]\n' > moderate.txt
zgrep -E "#|LOW" WC_outgroups.singl.outmiss.vcf.gz | bcftools query -f '%CHROM %POS  %REF  %ALT [\t%GT:%DP]\n' > low.txt

# Extract different profiles from filtered shared positions
zgrep -E "#|HIGH" WC_outgroup.singl.outmiss.intersect.vcf.gz | bcftools query -f '%CHROM %POS  %REF  %ALT [\t%GT:%DP]\n' > high.intersect.txt
zgrep -E "#|MODERATE" WC_outgroup.singl.outmiss.intersect.vcf.gz | bcftools query -f '%CHROM %POS  %REF  %ALT [\t%GT:%DP]\n' > moderate.intersect.txt
zgrep -E "#|LOW" WC_outgroup.singl.outmiss.intersect.vcf.gz | bcftools query -f '%CHROM %POS  %REF  %ALT [\t%GT:%DP]\n' > low.intersect.txt



