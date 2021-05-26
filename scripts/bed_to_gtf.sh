# convert bed file to gtf file for use with QTLtools quan
# also removes peaks on chr M, X, and Y
bed=$1
gtf=$2
sed '1d' $bed | \
awk 'BEGIN{OFS="\t"}{print $1, ".\texon",$2,$3,".\t.\t.\tgene_id\t\"peak" NR"\";"}' | \
 grep -v M | grep -v Y | grep -v X > $gtf
