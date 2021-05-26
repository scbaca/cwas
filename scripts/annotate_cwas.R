# annotate CWAS results with TWAS results and epigenomic data from LNCaP
library(dplyr)
library(gplots)

args <- commandArgs(trailingOnly = TRUE)
cwas.tested.file = args[1]
cwas.sig.file = args[2]
twas.file = args[3]
summary.file = args[4]
genes.file = args[5]
hichip.file = args[6]
lncap.ar.file = args[7]
lncap.H3K27ac.file = args[8]

twas=read.table(twas.file,header=T) #this file has only the models with cv pval < 0.05
#remove the few cases where there was a cross-validated model but no prediction:
twas=twas[!is.na(twas$TWAS.P),]
twas.sig=4.58e-7

cwas=read.table(cwas.tested.file,header=T) 
cwas.sig=read.table(cwas.sig.file,header=T)

#loop through significant cwas models and summarize info
win = 1000000

genes=read.table(genes.file, header=T, as.is=T)
genes$STRAND=ifelse(genes$STRAND=="+",1,0)

hc=read.table(hichip.file, header=T, as.is=T, skip=1)
colnames(hc)=c("CHR", "A0", "A1", "SCORE")

ar=read.table(lncap.ar.file, header=F, as.is=T)
colnames(ar)=c("CHR", "P0", "P1")

k27=read.table(lncap.H3K27ac.file, header=F, as.is=T)
colnames(k27)=c("CHR", "P0", "P1")

cwas.sig$TWAS.GENES = rep("none", nrow(cwas.sig))
cwas.sig$TWAS.GENES.HI = rep("none", nrow(cwas.sig))
cwas.sig$GENES = rep("none", nrow(cwas.sig))
cwas.sig$PROMOTER = rep("none", nrow(cwas.sig))
cwas.sig$HICHIP = rep("none", nrow(cwas.sig))
cwas.sig$LNCAP.AR = rep("FALSE", nrow(cwas.sig))
cwas.sig$LNCAP.H3K27ac = rep("FALSE", nrow(cwas.sig))

for(i in 1:nrow(cwas.sig)) {
	#nearby twas genes (all tissues):
	hits = twas$CHR==cwas.sig$CHR[i] & 
	  !((twas$P0<cwas.sig$P0[i]-win & twas$P1<cwas.sig$P0[i]-win) | (twas$P0>cwas.sig$P0[i]+win & twas$P1 > cwas.sig$P1[i]+win)) &
	  (twas$TWAS.P < twas.sig)
	if(any(hits)) cwas.sig$TWAS.GENES[i] = paste(unique(twas$ID[hits]), collapse="|")

	#nearby twas genes that have Z scores > 90% of the best GWAS Z score:
	hits = hits & (twas$TWAS.Z^2 > 0.9 * twas$BEST.GWAS.Z^2)
	if(any(hits)) cwas.sig$TWAS.GENES.HI[i] = paste(unique(twas$ID[hits]), collapse="|")

	#nearby genes:
	hits = genes$CHR==cwas.sig$CHR[i] &
	  (genes$TSS > cwas.sig$P0[i] - win & genes$TSS < cwas.sig$P1[i] + win)

	dist = cbind(cwas.sig$P0[i] - genes$TSS[hits], cwas.sig$P1[i] - genes$TSS[hits], genes$STRAND[hits]) %>% apply(1,FUN=function(x){
		ix = min(abs(x[1:2]))==abs(x[1]); #find closer peak boundary
		d = x[ifelse(ix,1,2)]; #find distance to closer peak boundary (retain sign)
		ifelse((x[3]==1 & d < 0) | (x[3]==0 & d > 0), -1*abs(d),abs(d))}) #adjust sign based on strand

	glist = paste(genes$GENE[hits], dist, sep=":")


	prom = glist[dist > -5000 & dist < 1000]
	if(length(prom) > 0) cwas.sig$PROMOTER[i] = paste(prom, collapse = "|")

	glist = glist[order(abs(dist))]

	if(any(hits)) cwas.sig$GENES[i] = paste(glist, collapse= "|")
	
	#HiChIP loops:
	hits0 = (hc$CHR==cwas.sig$CHR[i]) & (hc$A0 > cwas.sig$P0[i] - 5000 & hc$A0 < cwas.sig$P1[i] + 5000)
	hits1 = (hc$CHR==cwas.sig$CHR[i]) & (hc$A1 > cwas.sig$P0[i] - 5000 & hc$A1 < cwas.sig$P1[i] + 5000)

	hc.genes=character()
	for(j in which(hits0)) {
		g = genes$GENE[(genes$CHR==hc$CHR[j]) & (genes$TSS > hc$A1[j] - 5000 & genes$TSS < hc$A1[j] + 5000)]
		hc.genes=c(hc.genes, g)	
	}
	for(j in which(hits1)) {
		g = genes$GENE[(genes$CHR==hc$CHR[j]) & (genes$TSS > hc$A0[j] - 5000 & genes$TSS < hc$A0[j] + 5000)]
		hc.genes=c(hc.genes, g)	
	}	
	if(length(hc.genes)>0) cwas.sig$HICHIP[i]=paste(unique(hc.genes), collapse="|")

	#LNCaP AR overlap:
	hits = ar$CHR==cwas.sig$CHR[i] &
          !((ar$P0<cwas.sig$P0[i] & ar$P1<cwas.sig$P0[i]) | (ar$P0>cwas.sig$P0[i] & ar$P1 > cwas.sig$P1[i]))	
	if(any(hits)) cwas.sig$LNCAP.AR[i] = "TRUE"

        #LNCaP H3K27ac overlap:
        hits = k27$CHR==cwas.sig$CHR[i] &
          !((k27$P0<cwas.sig$P0[i] & k27$P1<cwas.sig$P0[i]) | (k27$P0>cwas.sig$P0[i] & k27$P1 > cwas.sig$P1[i]))
        if(any(hits)) cwas.sig$LNCAP.H3K27ac[i] = "TRUE"

}

write.table(cwas.sig,file=summary.file,sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
