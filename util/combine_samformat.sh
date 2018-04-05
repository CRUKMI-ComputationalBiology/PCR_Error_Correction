
infile=$1
num=$2
outfile=$3
statfile=$4

dum="_"

for((a=0; a<num; a++)) 
do
	cat $infile$dum$a.fa >> $outfile.sam
	rm $infile$dum$a.fa
done

samtools view -bS -o $outfile.bam $outfile.sam
#samtools view -bT $reference $outfile.sam > $outfile.bam

rm $outfile.sam
#mv "Stats_BaseCalls.txt" $outfile.stats
