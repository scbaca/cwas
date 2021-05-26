# reformat bam fields to work wtih ASEReadCounter
# get info from read names to add to the bams using picard AddOrReplaceReadGroups:
# this extracts the following and makes sure it's consistent for every read:
# instrument name : run id : flow cell id : lane
# also remove "chr" from sequence names and rmove chrM
bam_in=$1
bam_out=$2
name=$3
path=$4
dict=$5

mkdir -p $path

# get info from read names to add to the bams using picard AddOrReplaceReadGroups:
# (instrument name : run id : flow cell id : lane)
printf "adding read group info for %s\n" $bam_in

# Due to the formating of the NKI.Stello sample fastq files, the read group info is not present in the bams
# In these cases, will use dummy readgroup info

if [ $(samtools view $bam_in | head -n 1 | grep "^SRR" | wc -l) = 1 ]; then 
echo readgroup info missing for $in_bam. Adding placeholders.

	RID=RID:placeholder
	PU=PU:placeholder
else
	RID=$(samtools view $bam_in | head -n 1 | awk 'BEGIN{FS=":"}{print $1 ":" $2}')
	PU=$(samtools view $bam_in | head -n 1 | awk 'BEGIN{FS=":"}{print $3 ":" $4}')

fi

picard AddOrReplaceReadGroups INPUT=$bam_in OUTPUT=${path}/${name}.rg.tmp.bam RGID=$RID RGLB=lib1 RGPL=ILLUMINA RGPU=$PU:$name RGSM=$name

printf "removing chrM from $bam_in and deleting from chromosome names \n" $bam_in
samtools view -h ${path}/${name}.rg.tmp.bam | \
grep -v chrM | sed 's/chr//g' | samtools view -S -b > ${path}/${name}.nochr.tmp.bam

printf "reordering %s\n" $bam_in
picard ReorderSam I=${path}/${name}.nochr.tmp.bam O=$bam_out SD=$dict
samtools index $bam_out
rm ${path}/${name}.rg.tmp.bam ${path}/${name}.nochr.tmp.bam
