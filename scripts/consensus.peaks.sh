# Create a consensus peak set from all peak bed files

beds=$1
genome=$2
n=$3

sort -k1,1 -k2,2n $genome > analysis/consensus_peaks/genome.sorted

# partition genome in 50bp
bedtools makewindows -g analysis/consensus_peaks/genome.sorted -w 50 > analysis/consensus_peaks/win.50bp.bed
 
bedtools intersect -c -sorted -a analysis/consensus_peaks/win.50bp.bed -b $beds > analysis/consensus_peaks/counts.bed

#keep intervals with peaks in more than specified number of samples, buffer by 100bp, and merge:

awk -v n=$n 'BEGIN{OFS="\t"} $4>=n {print $1,$2,$3}' analysis/consensus_peaks/counts.bed | bedtools slop -i - -b 100 -g analysis/consensus_peaks/genome.sorted | bedtools merge -i - > analysis/consensus_peaks/consensus.bed

awk 'BEGIN{OFS="\t"; print "CHR\tP0\tP1\tNAME\tCENTER"}{print $1, $2, $3, $1 ":" $2 "-" $3, ($3-$2)/2+$2}' analysis/consensus_peaks/consensus.bed | sed 's/chr//' > analysis/consensus_peaks/consensus.peaks.stratas.tsv
