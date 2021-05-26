# perform FDR correction for stratAS results and identify peaks with significantly imbalanced SNPs. 
# also generates a number of qc/santiy check plots, most of which are not needed

args <- commandArgs( trailingOnly = TRUE )

library(qvalue)
library(ggplot2)
library(tidyr)
library(Hmisc)
library(ggplot2)

set.seed(123)

SIG=0.05 # Q value significance threshold
READ_THRESH=20 # SNPs must have at least this many reads used in AI calculation to be tested for significance
BAL_THRESH=100 # SNPs must have at least this many reads used in AI calculation to be added to list of covered but balanced SNPs
dir.create("analysis/summary/bed", recursive=T)

results.file=args[1]
samples.file=args[2]
plot.file=args[3]

cat("reading in stratAS results file", results.file, "\n")
all.results=read.table(results.file,sep="\t",header=TRUE,colClasses=c("character", "numeric", "character", "numeric", "numeric", "character", rep("numeric",10), rep("character",6)))

results=all.results[all.results$N.READS >= READ_THRESH,]
results$POSm1=results$POS-1

# get qvalues values
results$BBINOM.Q=qvalue(results$ALL.BBINOM.P)$qvalue
results$BBINOM.C0.Q=qvalue(results$C0.BBINOM.P)$qvalue
results$BBINOM.C1.Q=qvalue(results$C1.BBINOM.P)$qvalue
results$DIFF.Q=qvalue(results$DIFF.BBINOM.P)$qvalue

# find snps with significant AI in C0
sig.C0=subset(results,BBINOM.C0.Q<SIG)
sig.C0=sig.C0[order(sig.C0$BBINOM.C0.Q),]
write.table(sig.C0, "analysis/summary/sig.C0.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(unique(sig.C0[,c("CHR","P0","P1","NAME")]), "analysis/summary/bed/sig.C0.peak.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sig.C0[,c("CHR","POSm1","POS","RSID")], "analysis/summary/bed/sig.C0.snp.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# find snps with significant AI in C1
sig.C1=subset(results,BBINOM.C1.Q<SIG)
sig.C1=sig.C1[order(sig.C1$BBINOM.C1.Q),]
write.table(sig.C1, "analysis/summary/sig.C1.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(unique(sig.C1[,c("CHR","P0","P1","NAME")]), "analysis/summary/bed/sig.C1.peak.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sig.C1[,c("CHR","POSm1","POS","RSID")], "analysis/summary/bed/sig.C1.snp.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# find snps with significant AI in both conditions
sig.BOTH=subset(results,BBINOM.Q<SIG)
sig.BOTH=sig.BOTH[order(sig.BOTH$BBINOM.Q),]
write.table(sig.BOTH, "analysis/summary/sig.BOTH.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(unique(sig.BOTH[,c("CHR","P0","P1","NAME")]), "analysis/summary/bed/sig.BOTH.peak.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sig.BOTH[,c("CHR","POSm1","POS","RSID")], "analysis/summary/bed/sig.BOTH.snp.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# find snps with significantly differential AI across conditions
sig.DIFF=subset(results,DIFF.Q<SIG)
if (nrow(sig.DIFF) > 0) {
	sig.DIFF=sig.DIFF[order(sig.DIFF$DIFF.Q),]

	sig.DIFF$DELTA=sig.DIFF$C1.AF - sig.DIFF$C0.AF #difference in mu across states

	write.table(sig.DIFF, "analysis/summary/sig.DIFF.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
	write.table(unique(sig.DIFF[,c("CHR","P0","P1","NAME")]), "analysis/summary/bed/sig.DIFF.peak.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	write.table(sig.DIFF[,c("CHR","POSm1","POS","DELTA")], "analysis/summary/bed/sig.DIFF.snp.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# find a set of comparison SNPs with no imbalance wtih coverage > threshold
notsig=subset(results,!is.na(BBINOM.C0.Q) & !is.na(BBINOM.C1.Q) & !is.na(BBINOM.Q) & !is.na(DIFF.Q))
notsig=subset(notsig,BBINOM.C0.Q >= SIG & BBINOM.C1.Q >= SIG & BBINOM.Q >= SIG & DIFF.Q >= SIG & N.READS >= BAL_THRESH)
write.table(notsig, "analysis/summary/notsig.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

nr=hist(results$N.READS,breaks=c(1:max(results$N.READS)))
nh=hist(results$N.HET,breaks=c(1:max(results$N.HET)))

# for finding no. significant snps at given q value threshold:
C0.hist=hist(sig.C0$BBINOM.C0.Q,breaks=50)
C1.hist=hist(sig.C1$BBINOM.C1.Q,breaks=50)
BOTH.hist=hist(sig.BOTH$BBINOM.Q,breaks=50)

# for finding no. signfiicant peaks at given q value threshold:
peaks.C0=aggregate(sig.C0, by=list(sig.C0$NAME), FUN=min)
peaks.C1=aggregate(sig.C1, by=list(sig.C1$NAME), FUN=min)
peaks.BOTH=aggregate(sig.BOTH, by=list(sig.BOTH$NAME), FUN=min)

peaks.C0.hist=hist(peaks.C0$BBINOM.C0.Q,breaks=50)
peaks.C1.hist=hist(peaks.C1$BBINOM.C1.Q,breaks=50)
peaks.BOTH.hist=hist(peaks.BOTH$BBINOM.Q,breaks=50)

# plots
pdf("analysis/summary/stratas_plot.pdf")
plot(y=nr$counts,x=c(1:(max(results$N.READS)-1)),ylab="N SNPs",xlab="N reads used for AI calculation")
plot(y=nr$counts,x=c(1:(max(results$N.READS)-1)),ylab="N SNPs",xlab="N reads used for AI calculation", xlim=c(0,200))
plot(y=nh$counts,x=c(1:(max(results$N.HET)-1)),ylab="N SNPs",xlab="N heterozygous individuals used for AI calculation")
plot(y=nrow(results)-cumsum(nh$counts),x=nh$mids,ylab="N SNPs",xlab="N heterozygous individuals used for AI calculation (cumulative)")

plot(y=cumsum(C0.hist$counts),x=C0.hist$mids,ylab="N significant SNPS", xlab="Q-value cutoff", main="No. significant SNPs: C0")
plot(y=cumsum(C1.hist$counts),x=C1.hist$mids,ylab="N significant SNPS", xlab="Q-value cutoff", main="No. significant SNPs: C1")
plot(y=cumsum(BOTH.hist$counts),x=BOTH.hist$mids,ylab="N significant SNPS", xlab="Q-value cutoff", main="No. significant SNPs: Both conditions")

plot(y=cumsum(peaks.C0.hist$counts),x=peaks.C0.hist$mids,ylab="N significant peaks", xlab="Q-value cutoff", main="No. significant peaks: C0")
plot(y=cumsum(peaks.C1.hist$counts),x=peaks.C1.hist$mids,ylab="N significant peaks", xlab="Q-value cutoff", main="No. significant peaks: C1")
plot(y=cumsum(peaks.BOTH.hist$counts),x=peaks.BOTH.hist$mids,ylab="N significant peaks", xlab="Q-value cutoff", main="No. significant peaks: Both conditions")

n=nrow(results)
ix=sample(c(1:n),size=5000,replace=F)

plot(y=-1*log(results$BBINOM.Q[ix],10), x=results$N.READS[ix],ylab = "-log10(Q)", xlab = "N reads", main="AI significance vs coverage: combined conditions", ylim=c(0,20))
abline(h= -1*log(SIG,10),col="red")
plot(y=-1*log(results$BBINOM.Q[ix],10), x=results$N.READS[ix],ylab = "-log10(Q)", xlab = "N reads", main="AI significance vs coverage: combined conditions", xlim=c(0,200), ylim=c(0,20))
abline(h= -1*log(SIG,10),col="red")
plot(y=-1*log(results$DIFF.Q[ix],10), x=results$N.READS[ix],ylab = "-log10(Q)", xlab = "N reads", main="Differential AI significance vs coverage")

plot(y=sig.C0$C0.AF, x=sig.C0$C1.AF, ylab="mu in C0", xlab="mu in C1", main=paste0("mu in C0 vs C1 for snps with AI in C0 (q<", SIG,")"))
plot(y=sig.C1$C0.AF, x=sig.C1$C1.AF, ylab="mu in C0", xlab="mu in C1", main=paste0("mu in C0 vs C1 for snps with AI in C1 (q<", SIG,")"))

dev.off()

#significance vs read depth saturation plot:
mat=numeric() 
step=30
bins=seq(0,800,step)
for (i in bins) {
	ix=all.results$N.READS > i-step & all.results$N.READS < i
	n = sum(ix)
	if (n<200) next()
#	s = sum(results$BBINOM.Q[ix] < 0.05) / n
	s = qvalue(all.results$ALL.BBINOM.P[ix])$pi0
	mat <- rbind(mat, c(s,n,i))
}
mat <- as.data.frame(mat)
colnames(mat) = c("sig.fraction", "n.snps", "reads.bin")

ggplot(mat,aes(y=sig.fraction,x=reads.bin,size=n.snps)) + geom_point() + theme_classic() + scale_size(range = c(0, 5)) + ylim(0,1) + ggtitle(paste0("Estimated proportion of true null pvalues (pi0)", SIG, " vs read depth\n (only includes bins with >= 200 SNPs)"))
ggsave("analysis/summary/sig.vs.readdepth.pdf")

#test read depth cutoffs for significance testing on chr1
cuts=numeric() 
m=apply(cbind(all.results$N.READS, all.results$N.HET, all.results$ALL.BBINOM.P, all.results$CHR),1:2,as.numeric)
m = m[m[,4]==1,] #limit to chr1
for (r in seq(0,200,5)) {
	for (i in seq(1,10,1)) {
		cm = m[m[,1] >= r & m[,2] >= i,]
		c = qvalue(cm[,3])$qvalue
		nt = length(c)
		cuts <- rbind(cuts, c(r, i, nt, sum(c<SIG)))
	}
}

cuts <- as.data.frame(cuts)
colnames(cuts)=c("read.cutoff","indiv.cutoff","n.tested","n.sig")

ggplot(cuts,aes(y=n.sig,x=read.cutoff,size=n.tested,color=indiv.cutoff)) + geom_jitter() + theme_classic() + scale_size(range = c(0, 4)) + ggtitle(paste0("Number of significantly imbalanced SNPs at various read depth cutoffs for testing\n Data for chr1 only")) + scale_color_gradient(low="blue", high="gold")
ggsave("analysis/summary/cutoffs.analysis.pdf")

#get and plot individual level read counts
get_counts <- function(col) {
	lapply(col,function(x) as.numeric(unlist(strsplit(x,","))))
	}

IND.C0=get_counts(results$IND.C0)
IND.C0.COUNT.REF=get_counts(results$IND.C0.COUNT.REF)
IND.C0.COUNT.ALT=get_counts(results$IND.C0.COUNT.ALT)
IND.C1=get_counts(results$IND.C1)
IND.C1.COUNT.REF=get_counts(results$IND.C1.COUNT.REF)
IND.C1.COUNT.ALT=get_counts(results$IND.C1.COUNT.ALT)

#plot the distribution of coverage for snps (without combining multiple SNPs in a haplotype)
pdf("analysis/summary/snp.coverage.distribution.pdf")
snp.counts=c(unlist(IND.C0.COUNT.REF)+unlist(IND.C0.COUNT.ALT),unlist(IND.C1.COUNT.REF) + unlist(IND.C1.COUNT.ALT))
snp.counts.hist=hist(snp.counts[snp.counts<50],breaks=50,main="histogram of SNP read coverage", xlab="SNP read coverage") 

#plot AF for blaanced and imbalanced SNPs:
sig.hist=hist(results$ALL.AF[results$BBINOM.Q<SIG],breaks=40,main="histogram of AF for signficantly imbalanced SNPs",xlab="Allele fraction")

notsig.hist=hist(results$ALL.AF[results$BBINOM.Q>SIG & results$N.READS>BAL_THRESH],breaks=40,xlim=c(0,1),main="histogram of AF for balanced SNPs with high coverage",
xlab="Allele fraction")

dev.off()

#plot individual level read counts for selected snps
plotlist=list()
plotcount=1
plot.file="analysis/summary/allele_counts_plot.pdf"

names=read.table(samples.file, header=T, sep="\t", colClasses=(c("character","NULL")))$ID

#for now, just plotting most significant results:
snps.to.plot=head(order(results$BBINOM.Q),n=40)

for(i in snps.to.plot) {
	n0=length(unlist(IND.C0[i]))
	n1=length(unlist(IND.C1[i]))
	if(n0>0) {
		df0<-data.frame("rsid"=results$RSID[i], 
		  "condition"=rep("condition.0",n0),"ind"=unlist(IND.C0[i]),
		  "ref"=unlist(IND.C0.COUNT.REF[i]),"alt"=unlist(IND.C0.COUNT.ALT[i]))
	}
	if(n1>0) {
		df1<-data.frame("rsid"=results$RSID[i], 
		  "condition"=rep("condition.1",n1),"ind"=unlist(IND.C1[i]),
		  "ref"=unlist(IND.C1.COUNT.REF[i]),"alt"=unlist(IND.C1.COUNT.ALT[i]))
	}
	if(n0>0 && n1>0) { 
		plotmat <- rbind(df0,df1)
	} else if(n0>0 && n1==0) {
		plotmat <- df0
	} else {
		plotmat <- df1
	}
	plotmat <- plotmat %>% gather("ref","alt",key="allele",value="counts")

	plotmat$ind=factor(unlist(lapply(plotmat$ind,function(x){names[x]})), levels=names[order(names)])

	p <- ggplot(plotmat) + geom_bar(aes(x=allele,y=counts,fill=allele),stat="identity", position="stack", width=0.8, color="black")

	p <- p + facet_grid(~ind, drop=FALSE) +
	  theme_classic() +
	  labs(title=paste0("SNP ", results$RSID[i],
	  " chr", results$CHR[i], ":", results$POS[i] ,", region ", 
	  results$NAME[i]), subtitle=paste0("Condition 0 AF=", results$C0.AF[i],
	  ", p=",results$C0.BBINOM.P[i] ,"; Condition 1 AF=", results$C1.AF[i], 
	  ", p=",results$C1.BBINOM.P[i] ,"\np.diff=", results$DIFF.BBINOM.P[i] )) +
	  theme(axis.text.x = element_text(angle = 90)) +
	  theme(strip.text.x = element_text(angle = 90)) 		  
	plotlist[[plotcount]] <- p	
	plotcount <- plotcount + 1
}
pdf(plot.file, width=30,height=5)
for(i in seq_along(plotlist)) {
	print(plotlist[[i]])
}
dev.off()
