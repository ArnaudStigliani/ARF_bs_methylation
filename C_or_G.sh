#!/bin/bash


> GC_rate
> temp3

for sample in "ARF5_pos" "ARF5_neg" "ARF2_pos" "ARF2_neg"
do

> temp2
for position in "T1" "G2" "T3" "C4" "G5" "G6"
do
 
awk 'BEGIN{OFS="\t"}{print $5,$6,$7,$8}' ${position}_C_all_met_${sample}.bed | sed 's/-1$/1\t1\t-/' | sed 's/1$/1\t1\t+/'| bedtools getfasta -fi ~/Data/tair10/tair10.fas -fo temp.fas -bed -  -s
paste <(grep G temp.fas | wc -l) <(grep C temp.fas | wc -l) >>temp2
done
cat <(echo $sample"	"$sample) temp2 > temp3
paste  GC_rate temp3   > temp2
cat temp2 > GC_rate
done
rm temp.fas temp2 temp3

cat GC_rate

exit 0
