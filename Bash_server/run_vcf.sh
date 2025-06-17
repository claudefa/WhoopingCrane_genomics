while read chrom;
do
	sbatch snpAD_historical.sh $chrom
	sbatch snpAD_modern10.sh $chrom
	sbatch snpAD_modern5_15.sh $chrom
done < BigScaffolds_rohan.txt

while read line;
do
	chrom=$(echo $line | awk '{print $1}' )
	sample=$(echo $line | awk '{print $2}' )
       sbatch snpAD_modern10.sh $chrom $sample
done < snpad_SamplesChr.txt
