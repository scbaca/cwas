# convert gwas summary statistics files with rsid to bed file with hg19 snp positions
# includes only snps with abs(Z)>5.45
# of, if there are no SNPs meeting this criteria, take the 200 SNPs with the highest Z score

args=commandArgs(trailingOnly = TRUE)

snps.file=args[1]
bed.file=args[2]

t=read.table(text = gsub(":", " ", readLines("LDREF/hm3.pos")))
snps=read.table(snps.file,header=TRUE)

#limit to significant SNPs or top SNPs by Z-score:
if(any(abs(snps$Z)>5.45)) {
	snps=subset(snps, abs(Z)>5.45)
} else {
	snps=snps[order(-abs(snps$Z))[1:200],]
}

m = match(snps$SNP, t$V3)
bed = cbind(t[m,1],t[m,2]-1,t[m,2])
write.table(bed,bed.file,sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
