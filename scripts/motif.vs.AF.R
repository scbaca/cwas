#for het SNPs, plot difference in motif scores between the two alleles vs allelic imbalance

library(dplyr)
library(ggplot2)
library(cowplot)

args <- commandArgs( trailingOnly = TRUE )
motifs.by.allele=args[1]
ai.file=args[2]
motif.pdf=args[3]

m=read.table(motifs.by.allele,header=T,stringsAsFactors=F, sep = "\t")
ai=read.table(ai.file, sep = "\t", header=T)

m$CHR=m$chr
m$POS=m$pos

mer=merge(m,ai)
toplot = mer[,c("ref.score","alt.score","ALL.AF","motif")]

#set NAs to 0 for plotting (the NAs represent motifs that had such a low log odds score they were not detected by homer, (<0.001), so set them to 0)
toplot$alt.score[is.na(toplot$alt.score)]=0
toplot$ref.score[is.na(toplot$ref.score)]=0

toplot$diff=toplot$alt.score - toplot$ref.score
this.motif = unique(toplot$motif)
corr = cor.test (~ diff + ALL.AF, data = toplot, method="spearman")

ggplot(toplot, aes(y=diff,x=ALL.AF)) +
 geom_point (shape=16, alpha=0.3) + ylab("motif score difference (alt - ref)") +
 xlab("alt allele fraction") + geom_smooth(method='lm', se=F) +
 scale_x_continuous(breaks=c(0,0.5,1), limits = c(0,1)) +
 ggtitle(paste("difference in motif score vs AF \nfor ", nrow(toplot), " imbalanced SNPs\n",this.motif,"\np = ", corr$p.value, "\nrho = ", corr$estimate)) +
 theme_cowplot() +  theme(plot.title = element_text(size = 10))

ggsave(motif.pdf,width=3, height=3.5) 


