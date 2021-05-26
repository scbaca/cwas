#run ASEReadCounter to obtain allele-specific read counts

bam=$1
vcf=$2
genome=$3
out=$4

printf "running ASEReadCounter on reads in %s at snps in  %s\n" $bam $vcf 
gatk3 -T ASEReadCounter -I $bam \
                -sites $vcf \
                -o $out \
                -R $genome \
                -minDepth 0
