#!/bin/bash

#prepare stratAS input files

SAMPLE=$1
VCF=$2
COUNTS=$3
CNA=$4

STRATAS="stratAS"

mkdir -p analysis/stratas/${SAMPLE}

#get genotypes from vcf and create a sorted bed files with header: chr,pos0,pos1,rsid,genotype
bcftools query -f '%CHROM\t%POS0\t%END\t%ID\t%REF\t%FIRST_ALT[\t%GT]\n' $VCF | \
sort -k1,1 -k2,2n > analysis/stratas/${SAMPLE}/${SAMPLE}.snps.bed

#make a sorted bed file from the allele counts file with header chr,pos0,pos1,rsid,refcounts,altcounts
less $COUNTS | awk -v OFS='\t' '{print $1, $2-1, $2, $6, $7}' | \
sed '1d' | sort -k1,1 -k2,2n > analysis/stratas/${SAMPLE}/${SAMPLE}.counts.bed

echo merging genotypes and allele counts
bedtools intersect -wao -a analysis/stratas/${SAMPLE}/${SAMPLE}.snps.bed -b analysis/stratas/${SAMPLE}/${SAMPLE}.counts.bed > analysis/stratas/${SAMPLE}/${SAMPLE}.merged.bed

less analysis/stratas/${SAMPLE}/${SAMPLE}.merged.bed \
	| awk 'BEGIN {OFS="\t"; print "CHR\tPOS\tRSID\tREF\tALT\tHAP\tREF.READS\tALT.READS"} \
	{if($12==-1){print $1,$3,$4,$5,$6,$7,0,0} else {print $1,$3,$4,$5,$6,$7,$11,$12}}' \
	> analysis/stratas/${SAMPLE}/${SAMPLE}.counts

echo creating inp_counts file with covered hets only for stratas params.R 
less analysis/stratas/${SAMPLE}/${SAMPLE}.counts \
	| awk 'BEGIN {OFS=" "; print "CHR POS HAP REF.READS ALT.READS"} \
	NR>1 && $6 != "1|1" && $6 != "0|0" && $7+$8 > 0 {print $1, $2, $6, $7, $8}' \
	> analysis/stratas/${SAMPLE}/${SAMPLE}.counts.hets

echo preparing counts matrix file for stratAS
less analysis/stratas/${SAMPLE}/${SAMPLE}.counts | \
       awk 'BEGIN{OFS=" ";FS="\t"}{print $0}' | \
       sed 's/[|\t]/ /g' > analysis/stratas/${SAMPLE}/${SAMPLE}.counts.mat

# calculate paramaters (eg overdispersion paramaters) for stratas run
echo running params.R
Rscript $STRATAS/params.R \
	--min_snps 50 \
	--min_cov 5 \
	--inp_counts analysis/stratas/${SAMPLE}/${SAMPLE}.counts.hets \
	--inp_cnv ${CNA} \
	--group 10 \
	--out analysis/stratas/${SAMPLE}/${SAMPLE} \
	--id ${SAMPLE}
