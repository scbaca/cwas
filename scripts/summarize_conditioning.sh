#!/bin/bash
# summarize results from conditioning gwas significance on chromatin model SNPs

postprocess=$1
out=$2

cat $postprocess/*report | head -n1 | cut -f2-4,7-10 | sed 's/$/\tLOCUS/'> $out

prevchr=0
cat $postprocess/*report | grep -v FILE | sort -k2,2g -k3,3g | while read LINE; do

	#if this is the start of a new chromosome, reset the locus count
	chr=`echo $LINE | awk '{print $2}'`
	if [ $chr != $prevchr ]; then
		prevchr=$chr
		locus=1
	fi

	#repeat the entry the required number of times based on how many hits the locus contains
	n=`echo $LINE | awk '{print $5}'`
	for (( i = 1; i <= $n; i++ )); do
		echo $LINE | awk -v locus=$locus -v chr=$chr 'BEGIN{OFS="\t"}{print $2, $3, $4, $7, $8, $9, $10, "chr" chr ".loc_" locus}' >> $out
	done
	((locus++))
done

