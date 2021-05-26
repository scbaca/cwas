# Create a consensus peak set from all peak bed files:
# run this from the 'consensus' subfolder'

beds="../AR" #location of directory with bed files to merge
feat="AR" #name of the feature to merge

#for H3K27ac, use this:
beds="../H3K27ac" #location of directory with bed files to merge
feat="H3K27ac" #name of the feature to merge


echo using the following files:
ls $beds
cat ${beds}/*.bed | sort -k1,1 -k2,2n > ${feat}.all.bed

# partition genome in 50bp
if [ ! -f "hg19.genome" ]; then
	 mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo" | grep -v _ | grep -v size | sort -k1,1 -k2,2n > hg19.genome
fi

if [ ! -f "hg19.win.50bp.bed" ]; then
	bedtools makewindows -g hg19.genome -w 50 > hg19.win.50bp.bed
fi
 
if [ ! -f "${feat}.counts.bed" ]; then
	bedtools intersect -c -sorted -a hg19.win.50bp.bed -b ${feat}.all.bed > ${feat}.counts.bed
fi

#keep intervals with peaks in more than specified number of samples, buffer by 100bp, and merge:
n=`ls ${beds}/*.bed | wc -l | awk '{print $1 * 0.05}'`
echo keeping intervals with peaks in more than $n samples

awk -v n=$n 'BEGIN{OFS="\t"} $4>=n {print $1,$2,$3}' ${feat}.counts.bed | bedtools slop -i - -b 100 -g hg19.genome | bedtools merge -i - > ${feat}.consensus.bed

awk 'BEGIN{OFS="\t"; print "CHR\tP0\tP1\tNAME\tCENTER"}{print $1, $2, $3, $1 ":" $2 "-" $3, ($3-$2)/2+$2}' ${feat}.consensus.bed | sed 's/chr//' > ${feat}.consensus.stratas.bed
