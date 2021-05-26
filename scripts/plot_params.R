# plot the parameters used by stratas (eg, overdispersion parameters)
library(ggplot2)
library(cowplot)

args <- commandArgs( trailingOnly = TRUE )
loc.file=args[1]
glob.file=args[2]
loc.plot=args[3]
glob.plot=args[4]

loc=read.table(loc.file,header=T,colClasses=c("character",rep("NULL",3),rep("numeric",3),"NULL"))

p1 <- ggplot(loc,aes(y=CNV,x=ID)) + geom_jitter(pch=".") + geom_boxplot(outlier.shape=NA) + theme_cowplot() + ylab("CNV") + theme(axis.title.x=element_blank(), axis.text.x=element_blank()) 

p2 <- ggplot(loc,aes(y=MU,x=ID)) + geom_jitter(pch=".") + geom_boxplot(outlier.shape=NA) + theme_cowplot() + ylab("MU") + theme(axis.title.x=element_blank(), axis.text.x=element_blank()) + geom_hline(yintercept=0.5,color="blue")

p3 <- ggplot(loc,aes(y=PHI,x=ID)) + geom_jitter(pch=".") + geom_boxplot(outlier.shape=NA) + theme_cowplot() + ylab("PHI") + theme(axis.text.x= element_text(angle = 90, vjust=0.2, hjust=0.95))

plot_grid(p1,p2,p3,nrow=3, align="v")

ggsave(loc.plot, width=15, height=12)

glob=read.table(glob.file,header=T)
glob$fewSNPs = glob$N<1000

p1 = ggplot(glob,aes(y=MU,x=ID)) + geom_point() + theme_cowplot() + ylab("MU") + theme(axis.title.x=element_blank(), axis.text.x=element_blank()) + geom_hline(yintercept=0.5,color="blue")
p2 = ggplot(glob,aes(y=PHI,x=ID)) + geom_point() + theme_cowplot() + ylab("PHI") + theme(axis.title.x=element_blank(), axis.text.x=element_blank()) 
p3 = ggplot(glob,aes(y=N,x=ID,color=fewSNPs)) + geom_point() + theme_cowplot() + ylab("N SNPs") + ylim(0,50000) + theme(axis.text.x= element_text(angle = 90, vjust=0.2, hjust=0.95)) + scale_color_manual(values = c("black", "red")) + theme(legend.position="none")

plot_grid(p1,p2,p3,nrow=3, align="v")

ggsave(glob.plot, width=15, height=12)
