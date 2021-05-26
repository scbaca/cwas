#plot allele fraction (mu) vs qtl effect size (beta) for peaks with a cQTL and allelic imbalance
library(ggplot2)
library(Hmisc)
library(cowplot)

args <- commandArgs( trailingOnly = TRUE )
ai_file=args[1]
qtl_file=args[2]
plot=args[3]

ai=read.table(ai_file, sep="\t", header=T)
qtl=read.table(qtl_file, sep="\t", header=T, row.names=NULL)

ai$ID=paste(ai$CHR, ai$POS)
qtl$ID=paste(qtl$chr, qtl$snp.pos)

m=merge(ai,qtl,by="ID")

corr = cor.test (~ ALL.AF + beta, data = m, method="spearman")
ggplot(m, aes(x=ALL.AF,y=beta)) + geom_point(alpha=0.5, size=0.1, color="black") + geom_smooth(method='lm') + 
	ylab("qtl effect size (beta)") + 
	xlab("allele fraction (mu)") + 
	ggtitle(paste("AF vs beta for signifcantly imbalanced SNPs\np = ",corr$p.value, "\nrho = ", corr$estimate)) + 
	theme_cowplot()
ggsave(plot, height=3.5, width=3)
