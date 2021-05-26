# combine significane tests from qtl and ai for each snp using Stouffer's method
library(ggplot2)
library(Hmisc)
library(qvalue)
library(cowplot)

args <- commandArgs( trailingOnly = TRUE )
ai_file=args[1]
qtl_file=args[2]
out.comb=args[3]
out.snp=args[4]
out.bed=args[5]

ai=read.table(ai_file, sep="\t", header=T)
qtl=read.table(qtl_file, sep="\t", header=T, row.names=NULL)

ai$ID=paste(ai$CHR, ai$POS)
qtl$ID=paste(qtl$chr, qtl$snp.pos)

m=merge(ai,qtl,by="ID")

m$ai.z=qnorm(m$ALL.BBINOM.P/2,lower.tail=F)
ai.down=m$ALL.AF < 0.5
m$ai.z[ai.down] = -1 * m$ai.z[ai.down]

m$qtl.z=qnorm(m$pval/2,lower.tail=F)
qtl.down=m$beta < 0
m$qtl.z[qtl.down] = -1 * m$qtl.z[qtl.down]

m$comb.z = (m$ai.z + m$qtl.z) / sqrt(2)
m$comb.p = 2*(pnorm(abs(m$comb.z) , lower.tail=F))

q=qvalue(m$comb.p)
m$comb.q = q$qvalues

#write snps to file:
write.table(m, out.comb, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#output peaks and snps that are significant in combined test:
sig.comb=subset(m,comb.q<0.05)
sig.comb=sig.comb[order(sig.comb$comb.q),]
write.table(sig.comb, out.snp, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(unique(sig.comb[,c("CHR","P0","P1","NAME")]), out.bed, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
