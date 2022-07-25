# quantify overlap of foreground features (ie cQTL peaks) with targets (ie eQTL SNPs)

fg=$1
targets=$2
out=$3

fg_tot=`cat $fg | sort | uniq | awk '{ print $3-$2 }' | awk -f scripts/sum.awk` #sum basepairs in sig peaks
fg_hit=`cat $fg | sed 's/chr//' | sort -k1,1 -k2,2n | uniq | bedtools intersect -wa -a - -b $targets | sort | uniq | awk '{ print $3-$2 }' | awk -f scripts/sum.awk` # sum basepairs in sig peaks that overlap a risk snp
fg_num=`cat $fg | sort | uniq | wc -l` # num sig peaks

fg_enrich=`echo $fg_hit $fg_tot | awk '{ print 1e3 * $1/$2 }'` # portion of sig peak territory overlapping risk snp (x 1000)

printf "fg_tot\tfg_hit\tfg_num\tfg_enrich\n" > $out
printf "$fg_tot\t$fg_hit\t$fg_num\t$fg_enrich\n" >> $out


