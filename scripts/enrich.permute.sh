#!/bin/bash

n=$1 #batch number for this run
bg=$2 #background file (eg, bed file with all H3K27ac tested for imbalance)
fg=$3 #target (eg, imbalanced)  peaks
targets=$4 #bed file with eqtl snps
out_dir=$5

fg_num=`cat $fg | sort | uniq | wc -l` # num sig peaks

printf '' > ${out_dir}/group.${n}.bg.enrich

tmp=`echo ${out_dir}/bg.${n}.tmp`
for i in `seq 1 100`; do
	echo processing iteration $i for batch $n

        cat $bg | grep -v CHR | shuf -n $fg_num | sed 's/chr//' | sort -k1,1 -k2,2n > $tmp
	head -n 2 $tmp #tmp
        bg_tot=`cat $tmp | sort | uniq | awk '{ print $3-$2 }' | awk -f scripts/sum.awk`
        bg_hit=`bedtools intersect -wa -a $tmp -b $targets | sort | uniq | awk '{ print $3-$2 }' | awk -f scripts/sum.awk`
        bg_enrich=`echo $bg_hit $bg_tot | awk '{ print 1e3 * $1/$2 }'` 

	echo $bg_enrich >> ${out_dir}/group.${n}.bg.enrich
done
rm $tmp
