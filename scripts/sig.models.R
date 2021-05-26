# take fusion output run on all models (including those with non-significant cross-validation p-values) and apply filtering steps and multiple hypothesis testing adjustment.
 
#Filtering approaches:
#lenient: take min p for each peak. keep if p < 0.05
#stringent: take min p for each peak. FDR correct for all models tested (ie, with non-NA results)
#gwas-centric: take all peaks with TWAS.P < 5e-8, restrict 
#gwas-restricted: only test peaks falling near top gwas snps

library(dplyr)
library(qvalue)
library(ggplot2)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
fusion.path=args[1]

dir.create(paste0(fusion.path,"/sig"))                                         
dir.create(paste0(fusion.path,"/cv"))                                         

files=list.files(fusion.path, pattern="fusion.*.txt", full.name=T)
m = do.call(rbind, lapply(files, function(x) read.table(x, stringsAsFactors=F,header=T, sep="\t",
  colClasses=c(rep("character",3), rep("numeric",4), "character", "numeric", "character", rep("numeric", 5), rep("character",3), rep("numeric",2)))))

# replace each MODELCV.PV with the minimum of the two values it contains (corresponding to allele-specific and qtl cross-validation performance)
p.tmp = matrix(as.numeric(unlist(str_split(m$MODELCV.PV, ","))), nrow=nrow(m), ncol=2, byrow=T)
p.best=suppressWarnings(apply(p.tmp, 1, FUN=function(x) min(x,na.rm=T))) #introduces Inf where both rows are NA
p.best[!is.finite(p.best)] = NA
m$MODELCV.PV=p.best

r2.tmp = matrix(as.numeric(unlist(str_split(m$MODELCV.R2, ","))), nrow=nrow(m), ncol=2, byrow=T)
m$MODELCV.R2=sapply(seq(1:nrow(r2.tmp)),FUN=function(x) ifelse(!is.na(p.tmp[x,1]) & p.best[x]==p.tmp[x,1], r2.tmp[x,1], r2.tmp[x,2]))

m = subset(m, !is.na(MODELCV.PV) & !is.na(BEST.GWAS.Z)) #remove models with no crossvalidation or no GWAS SNP

cat(paste0("total models considered: ", nrow(m), "\n"))

#lenient filtering (accepting nominally significant CV pvals):
m.lenient = subset(m,MODELCV.PV < 0.05)
m.lenient = m.lenient %>% group_by(ID) %>% slice(which.min(MODELCV.PV))
cat(paste0("total models with significant CV pval without fdr correction: ", nrow(m.lenient), "\n"))

m.lenient.sig = subset(m.lenient, TWAS.P < 0.05 / nrow(m.lenient))
cat(paste0("significant CWAS from these models after bf correction: ", nrow(m.lenient.sig), "\n\n"))

if(nrow(m.lenient) > 0) {
ggplot(m.lenient, aes(x=MODEL)) + geom_bar() + theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(paste0(fusion.path,"/cv/model.lenient.pdf"),heigh=3,width=3)
}

if(nrow(m.lenient.sig) > 0) {
ggplot(m.lenient.sig, aes(x=MODEL)) + geom_bar() + theme_classic() + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(paste0(fusion.path,"/sig/model.lenient.pdf"),heigh=3,width=3)
}

#more stringent filtering (select top model type, then fdr correct):
m.best = m %>% group_by(ID) %>% slice(which.min(MODELCV.PV))
cat(paste0("total models after selecting most significant model for a given peak: ", nrow(m.best), "\n"))
m.best = subset(m.best, qvalue(MODELCV.PV, pi0=1)$qvalues < 0.05)

cat(paste0("total models with significant CV pval using fdr correction after selecting top model type: ", nrow(m.best), "\n"))

m.best.sig = subset(m.best, TWAS.P < 0.05 / nrow(m.best))
cat(paste0("significant CWAS from these models after bf correction: ", nrow(m.best.sig), "\n\n"))

if(nrow(m.best) > 0) {
ggplot(m.best, aes(x=MODEL)) + geom_bar() + theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(paste0(fusion.path,"/cv/model.best.pdf"),heigh=3,width=3)
}

if(nrow(m.best.sig) > 0) {
ggplot(m.best.sig, aes(x=MODEL)) + geom_bar() + theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))                                     
ggsave(paste0(fusion.path,"/sig/model.best.pdf"),heigh=3,width=3)
}

#most stringent filtering (fdr correct all models tested):
m.stringent = subset(m,qvalue(MODELCV.PV, pi0=1 )$qvalue < 0.05)
m.stringent = m.stringent %>% group_by(ID) %>% slice(which.min(MODELCV.PV))
cat(paste0("total models with significant CV pval using fdr correction of all models considered: ", nrow(m.stringent), "\n"))

m.stringent.sig = subset(m.stringent, TWAS.P < 0.05 / nrow(m.stringent))
cat(paste0("significant CWAS from these models after bf correction: ", nrow(m.stringent.sig), "\n\n"))

if(nrow(m.stringent) > 0) {
ggplot(m.stringent, aes(x=MODEL)) + geom_bar() + theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(paste0(fusion.path,"/cv/model.stringent.pdf"),heigh=3,width=3)
}

if(nrow(m.stringent.sig) > 0) {
ggplot(m.stringent.sig, aes(x=MODEL)) + geom_bar() + theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  
ggsave(paste0(fusion.path,"/sig/model.stringent.pdf"),heigh=3,width=3)
}

#gwas-centric filtering (start with only signficant regions, keep only those with CV pvalues passing FDR correction):
m.gwas = subset(m, TWAS.P < 5e-8) 
if(nrow(m.gwas) > 0) {
	m.gwas = subset(m.gwas, qvalue(MODELCV.PV, pi0=1 )$qvalue < 0.05)
	m.gwas = m.gwas %>% group_by(ID) %>% slice(which.min(MODELCV.PV))
	cat(paste0("total models with significant CV pval using gwas-centric approach: ", nrow(m.gwas), "\n"))

}
m.gwas.sig = subset(m.gwas, TWAS.P < 0.05 / nrow(m.gwas))
cat(paste0("significant CWAS from these models after bf correction: ", nrow(m.gwas.sig), "\n\n"))


if(nrow(m.gwas.sig) > 0) {
	ggplot(m.gwas.sig, aes(x=MODEL)) + geom_bar() + theme_classic() +
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  
	ggsave(paste0(fusion.path,"/sig/model.gwas.pdf"))
}

#write files with siginificant CWAS peaks:
write.table(m.lenient.sig, paste0(fusion.path,"/sig/cwas.sig.lenient.txt"),sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(m.best.sig, paste0(fusion.path,"/sig/cwas.sig.best.txt"),sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(m.stringent.sig, paste0(fusion.path,"/sig/cwas.sig.stringent.txt"),sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(m.gwas.sig, paste0(fusion.path,"/sig/cwas.sig.gwas.txt"),sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#write files with all models with significant crossvalidation
write.table(m.lenient, paste0(fusion.path,"/cv/cwas.cv.lenient.txt"),sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE) 
write.table(m.best, paste0(fusion.path,"/cv/cwas.cv.best.txt"),sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(m.stringent, paste0(fusion.path,"/cv/cwas.cv.stringent.txt"),sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(m.gwas, paste0(fusion.path,"/cv/cwas.cv.gwas.txt"),sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#gwas-restricted filtering (only test peaks falling near top gwas regions)
gwas=str_replace(fusion.path,"/merged","/gwas.bed")

system(paste0("sed '1d' ", fusion.path, "/cv/cwas.cv.lenient.txt | cut -f4-6 | bedtools intersect -wa -a - -b ", gwas, " > ", fusion.path, "/cv/cwas.cv.gwas.restricted.txt"))

g=read.table(paste0(fusion.path, "/cv/cwas.cv.gwas.restricted.txt"))
incl=paste0("chr", g[,1], ":", g[,2], "-", g[,3])

m.gwas.restricted=subset(m.best, ID %in% incl)
cat(paste0("total models with nominally significant CV pval when limiting to gwas regions: ", nrow(m.gwas.restricted), "\n"))

m.gwas.restricted.sig = subset(m.gwas.restricted, TWAS.P < 0.05 / nrow(m.gwas.restricted))
cat(paste0("significant CWAS from these models after bf correction: ", nrow(m.gwas.restricted.sig), "\n\n"))

write.table(m.gwas.restricted.sig, paste0(fusion.path,"/sig/cwas.sig.gwas.restricted.txt"),sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
