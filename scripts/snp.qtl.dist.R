# plot distance from significant qtl peaks to their top cqtl snp 
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
file=args[1]

t=read.table(file)
colnames(t)=c("peak","start","end","chr", "dist","qval")
t$mid=(t$end-t$start)/2 + t$start
t$dist=t$dist+(t$start-t$mid)

t$inpeak = (t$mid + t$dist > t$start) & (t$mid + t$dist < t$end)
message("in peak: ", sum(t$inpeak)/nrow(t))

t$tenkbpeak = abs(t$dist) < 10000
message("within 10kb: ", sum(t$tenkbpeak)/nrow(t))

t$twentyfivekbpeak = abs(t$dist) < 25000
message("within 25kb: ", sum(t$twentyfivekbpeak)/nrow(t))

t$onehundredkbpeak = abs(t$dist) < 100000
message("within 100kb: ", sum(t$onehundred)/nrow(t))

t$twohundredbppeak = abs(t$dist) < 200
message("within 200bp: ", sum(t$twohundredbppeak)/nrow(t))

topquant = quantile(t$q, seq(0,1,0.05))[2]
message("peaks in within 25kb and in the top 5th percentile of q values: ", sum(t$q < topquant & t$twentyfivekbpeak)/sum(t$q < topquant))

message("wilcoxon p value for comparison of q values within and outside of 25kb: ", wilcox.test(x=t$q[t$twentyfivekbpeak], y=t$q[!t$twentyfivekbpeak])$p.value)

pdf("analysis/qtl/summary/snp.dist.pdf", width=3, height=3)
hist(t$dist,breaks=100, col="darkgray", border="darkgray", main="distance from peak center to top SNP", xlab="distance (bp)", ylab="number of SNPs", axes=F)
axis(side=1, at=seq(-1e6, 1e6, by=1e6))

hist(t$dist,breaks=100, col="darkgray", border="darkgray", main="distance from peak center to top SNP", xlab="distance (bp)", ylab="number of SNPs", axes=F)
axis(side=1, at=seq(-1e6, 1e6, by=1e6))
abline(v=c(-100000, 100000), col="blue", lty=2)

hist(t$dist,breaks=100, col="darkgray", border="darkgray", main="distance from peak center to top SNP", xlab="distance (bp)", ylab="number of SNPs", axes=F)
axis(side=1, at=seq(-1e6, 1e6, by=1e6))
abline(v=c(-25000, 25000), col="blue", lty=2)

hist(t$dist[t$dist > -100000 & t$dist < 100000],breaks=100, col="darkgray", border="darkgray", main="distance from peak center to top SNP", xlab="distance (bp)", ylab="number of SNPs", xlim=c(-100000,100000), axes=F)
axis(side=1, at=seq(-100000, 100000, by=100000))

hist(t$dist[t$dist > -1000 & t$dist < 1000],breaks=100, col="darkgray", border="darkgray", main="distance from peak center to top SNP", xlab="distance (bp)", ylab="number of SNPs", xlim=c(-1000,1000), axes=F)
abline(v=c(-200, 200), col="blue", lty=2)
axis(side=1, at=seq(-1000, 1000, by=1000))

plot(y= -log(t$qval,10), x=t$dist, main="qval vs distance to top snp", xlab="distance (bp)", ylab="-log10(qval)", pch=".")

plot(y= -log(t$qval,10), x=t$dist, main="qval vs distance to top snp", xlab="distance (bp)", ylab="-log10(qval)", pch=".")
abline(v=c(-25000, 25000), col="blue", lty=2)

plot(y= -log(t$qval,10), x=t$dist, main="qval vs distance to top snp", xlab="distance (bp)", ylab="-log10(qval)", pch=".")
abline(v=c(-100000, 100000), col="blue", lty=2)

dev.off() 
