#generate manhattan plots and zz plots for CWAS models
#also generate manhattan plot for TWAS models

library(qqman)
library(ggplot2)
library(dplyr)
library(reshape2)
library(gplots)
library(qvalue)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
fusion.path=args[1]
manhattan.file=args[2]
zzplot.base=args[3]

twas.file="twas_data/ONCOARRAY.TWAS.txt"

pdf(manhattan.file, height=3, width=6)

#add dummy snps at begining and end of chromosomes so manhattan plot will display with correct chr lengths:
d=readRDS("misc/dummysnps.rds")
#for (mod in c("lenient", "best", "stringent", "gwas")) {
for (mod in "best") {
	fusion.file=paste0(fusion.path,"/cv/cwas.cv.",mod,".txt")
	model=read.table(fusion.file,sep="\t",header=T, as.is=T)
	fusion.file.sig=paste0(fusion.path,"/sig/cwas.sig.",mod,".txt")
	model.sig=read.table(fusion.file.sig,sep="\t",header=T, as.is=T)
	sigline = max(model.sig$TWAS.P)
	sigline = ifelse(is.finite(sigline),sigline,0)

	plot=data.frame(SNP=model$BEST.GWAS.ID, CHR=model$CHR, P0=model$P0, P1=model$P1, BP=model$P0+(model$P1-model$P0)/2, P=as.numeric(model$TWAS.P), zscore=as.numeric(model$TWAS.Z), gwas.z = model$BEST.GWAS.Z, model=mod, file = model$FILE, modeltype=model$MODEL)
	plot$sig=plot$P <= sigline
	cat(paste0("significant snps for ", mod," model: ",sum(plot$sig),"\n"))
	cat(paste0("pval cutoff ", mod," model: ", sigline,"\n"))
	cat(paste0("zscore cutoff ", mod," model: ", min(abs(plot$zscore[plot$sig])), "\n\n"))

	plot=rbind(plot,d)
	if(!all(is.finite(log(plot$P,10)))) {ymax=200} else {ymax=max(-1*log(plot$P,10)) + 5}

	#plot all CWAS SNPs
        manhattan(plot, genomewideline = -log10(sigline), suggestiveline=F,                                        
                main=paste0(mod, " (", sum(plot$sig), " / ", nrow(plot), " peaks significant)"),                   
		ylim=c(0,ymax))

	#plot highly sigificant SNPs separately
	manhattan(plot, genomewideline = -log10(sigline), suggestiveline=F, 
		main=paste0(mod, " (", sum(plot$sig), " / ", nrow(plot), " peaks significant)"), 
		ylim=c(40,ymax))

	#plot less significant SNPs separately
        manhattan(plot, genomewideline = -log10(sigline), suggestiveline=F,                                   
                main=paste0(mod, " (", sum(plot$sig), " / ", nrow(plot), " peaks significant)"),              
                ylim=c(0,40)) 

	anysig=any(plot$sig)

	#remove dummy snps (that were used for manhattan plot)
	plot = subset(plot, SNP!="dummy")

	if(anysig) {
		colors=c(brewer.pal(3, name="Reds"), brewer.pal(3, name="Blues"))
		plot$modeltype = factor(plot$modeltype, 
		  levels = c("lasso", "lasso.as", "lasso.plasma", "top1.qtl", "top1.as", "top1"))
		names(colors) = levels(plot$modeltype)
		ggplot(plot, aes(x=abs(gwas.z), y=abs(zscore), col=modeltype)) + 
		  geom_point() + theme_classic() + 
		  ylim(2,25) + scale_color_manual(name="modeltype", values=colors) +
		  xlim(2,25) + ylab("abs(CWAS.Z)") + xlab("abs(best GWAS.Z)") + 
		  geom_vline(xintercept=5.45, linetype="dotted") +
		  geom_hline(yintercept=min(abs(plot$zscore[plot$sig])), linetype="dashed")
	}	

	ggsave(paste0(zzplot.base,".",mod,".pdf"), height=3, width=4)

	#how many peaks are significant by cwas without a significant gwas snp?
	message("number of peaks that are signfiicant by cwas without a signfiicant gwas snp:", sum(abs(plot$gwas.z) < 5.45 & plot$sig==T))

}

#plot twas data
twas=read.table(twas.file, header=T)
colnames(twas)[c(3,4,5,19,20)] = c("SNP", "CHR", "BP", "zscore", "P")
twas = subset(twas, !is.na(P)) #remove the few rows with NA p values
manhattan(twas, genomewideline = -log10(0.05/nrow(twas)), suggestiveline=F, main="TWAS")

dev.off()
