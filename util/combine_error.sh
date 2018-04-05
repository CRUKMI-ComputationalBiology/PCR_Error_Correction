
infile=$1
num=$2
outfile=$3

dum="_"

for((a=0; a<num; a++)) 
do
	cat $infile$dum$a >> $outfile
	rm $infile$dum$a
done


