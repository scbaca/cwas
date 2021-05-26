#get motif scores at the specified motif for sequence containing reference vs alternate SNP alleles

TSNPS=$1 #stratas output file containing SNP info for target (imbalanced) SNPS
COUNTS=$2 #counts matrix file to get ref/alt allele identities from - ie, all.counts.mat 
WIN=$3 #how many bp in each direction to expand BP search
OUT=$4 #output directory for homer results
GENOME=$5 #fasta file with ref genome sequence
MOTIFS=$6 #motifs to be tested

mkdir -p $OUT

echo getting allele info
awk 'BEGIN{OFS=" ";FS=" "} NR==FNR{a[$1,$2,$3]=$4 OFS $5; next} {FS="\t"} FNR>1 {print $1, $2, $3, a[$1,$2,$3]}' $COUNTS $TSNPS | sort -u -k1,4 > ${OUT}/alleles.tmp

# the sort uniq ensures that if there are more than two alleles at one locus, the one will be arbitrarily chosen. See eg rs116442826
awk -v w=$WIN 'BEGIN{OFS="\t";FS=" "}{print $1, $2-w-1, $2+w, $2, $4, $5, $3}'  ${OUT}/alleles.tmp > ${OUT}/snp.bed.win.tmp

echo creating sequences with ref allele
sed 's/chr//' ${OUT}/snp.bed.win.tmp | bedtools getfasta -fi $GENOME -bed - > ${OUT}/ref.fa

echo creating sequences with alt allele
awk -v w=$WIN 'NR==FNR{allele[NR]=$6; next} {if (FNR % 2 == 0) {start = substr($0,1,w) ; end=substr($0, w+2, w); print start allele[FNR/2] end} else {print $0}}' ${OUT}/snp.bed.win.tmp ${OUT}/ref.fa > ${OUT}/alt.fa

#find motifs, and remove any that don't overlap the SNP in the center of the window:
findMotifs.pl ${OUT}/ref.fa fasta ${OUT}/homer -find $MOTIFS | sort | tail -n +2 | \
awk '{if(($5=="+" && $2 > -1*length($3) && $2 <= 0) || \
($5=="-" && $2 >= 0 && $2 <= length($3))){print $0}}' > ${OUT}/ref.motifs

findMotifs.pl ${OUT}/alt.fa fasta ${OUT}/homer -find $MOTIFS | sort | tail -n +2 | \
awk '{if(($5=="+" && $2 > -1*length($3) && $2 <= 0) || \
($5=="-" && $2 >= 0 && $2 < length($3))){print $0}}' > ${OUT}/alt.motifs

#combine ref and alt motifs
cat ${OUT}/ref.motifs ${OUT}/alt.motifs | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$5}' | sort | uniq > ${OUT}/all.motifs

#make a temporary motif file with log odds threshold set low (to 0.0001):
sed 's/)\t.*/)\t0.0001/' $MOTIFS > ${OUT}/loose.motifs

#detect motifs in the ref and alt sequences using the "loose" motif thresholds
findMotifs.pl ${OUT}/ref.fa fasta ${OUT}/homer -find ${OUT}/loose.motifs | sort | tail -n +2 | \
awk '{if(($5=="+" && $2 > -1*length($3) && $2 <= 0) || \
($5=="-" && $2 >= 0 && $2 <= length($3))){print $0}}' > ${OUT}/ref.loose.motifs

findMotifs.pl ${OUT}/alt.fa fasta ${OUT}/homer -find ${OUT}/loose.motifs | sort | tail -n +2 | \
awk '{if(($5=="+" && $2 > -1*length($3) && $2 <= 0) || \
($5=="-" && $2 >= 0 && $2 < length($3))){print $0}}' > ${OUT}/alt.loose.motifs

awk 'BEGIN{OFS="\t"} NR==FNR{a[$1,$2,$4,$5]=$3 OFS $6; next} {if(length(a[$1,$2,$3,$4])>0){print $1, $2, $3, $4, a[$1,$2,$3,$4]} else {print $1, $2, $3, $4, "NA", "NA"}}' ${OUT}/ref.loose.motifs ${OUT}/all.motifs > ${OUT}/ref.motifs.final

awk 'BEGIN{OFS="\t"} NR==FNR{a[$1,$2,$4,$5]=$3 OFS $6; next} {if(length(a[$1,$2,$3,$4])>0){print $1, $2, $3, $4, a[$1,$2,$3,$4]} else {print $1, $2, $3, $4, "NA", "NA"}}' ${OUT}/alt.loose.motifs ${OUT}/all.motifs > ${OUT}/alt.motifs.final

#get SNP coordinates:
sed 's/:/\t/' ${OUT}/all.motifs | sed 's/-.*//' | awk -v w=$WIN 'BEGIN{OFS="\t"}{print $1, $2+w+1}' > ${OUT}/coords.tmp

printf "chr\tpos\tpeak\toffset\tmotif\tstrand\tref.seq\tref.score\talt.seq\talt.score\n" > ${OUT}/motifs.by.allele
paste ${OUT}/coords.tmp ${OUT}/all.motifs ${OUT}/ref.motifs.final ${OUT}/alt.motifs.final | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$11,$12,$17,$18}' >> ${OUT}/motifs.by.allele

rm ${OUT}/coords.tmp  ${OUT}/alleles.tmp ${OUT}/snp.bed.win.tmp
