# compare average Z-squared of cwas associations across different gwas studies
library(ggplot2)

args=commandArgs(trailingOnly = TRUE)
plotfile="analysis/fusion/compare_gwas/compare.gwas.pdf"
n=length(args)
names=gsub("analysis/fusion/","",args)
names=gsub("/merged/cv/cwas.cv.best.txt","",names)

sem = function(x) {sd(x)/sqrt(sum(!is.na(x)))}

zs = data.frame(names=names, z_squared=rep(0,n), sem=rep(0,n), nmodels=rep(0,n))
for (i in 1:n) {
	m=read.table(args[i],header=T)
	zs$z_squared[i] = mean(m$TWAS.Z^2)
	zs$sem[i] = sem(m$TWAS.Z^2)
	zs$nmodels[i] = nrow(m)
}

zs$z_squared=zs$z_squared/mean(zs$z_squared)
zs$sem=zs$sem/mean(zs$z_squared)

names=as.factor(names)
ggplot(zs, aes(y=z_squared,x=reorder(names, -z_squared), ymin=z_squared-sem, ymax=z_squared+sem)) + geom_col() + geom_linerange() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + xlab("") + ylab("mean z-squared")

ggsave(plotfile, height=5, width=4)
