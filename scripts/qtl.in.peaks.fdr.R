#get q-values for qtl SNPs in peaks
library(qvalue)

args <- commandArgs( trailingOnly = TRUE )
qtl_file=args[1]
out.snp=args[2]
out.bed=args[3]

qtl=read.table(qtl_file, sep="\t", header=T, row.names=NULL)
q=qvalue(qtl$pval)
qtl$q=q$qvalues

qtl = subset(qtl, q < 0.05)

#output peaks and snps that are significant by q-val:
write.table(qtl, out.snp, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(unique(qtl[,c("chr","start","end")]), out.bed, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
