# Use copywriteR package to estimate copy number profiles form ChIP-seq data using off-target reads

library(CopywriteR)
library(matrixStats)
library(gtools)
library(data.table)
library(chipseq)
library(GenomicAlignments)
library(Rsamtools)
library(futile.logger)
library(CopywriteR)
library(IRanges)
library(S4Vectors)
library(DNAcopy)
library(GenomicRanges)
library(GenomeInfoDb)
library(futile.logger)

source("scripts/CopywriteRCustomBed.R")
source("scripts/private.R")

args <- commandArgs( trailingOnly = TRUE )

bam.path=args[1]
name=args[2]
cna.path=args[3]
bed=args[4]

ncores=8

preCopywriteR(output.folder = file.path(cna.path), bin.size = 400000, ref.genome = "hg19")

load(file = file.path(cna.path, "hg19_400kb", "blacklist.rda"))

bp.param <- MulticoreParam(workers=ncores)

samples <- list.files(path = bam.path, pattern = ".bam$", full.names = TRUE)
controls <- samples
sample.control <- data.frame(samples, controls)

CopywriteRCustomBed(sample.control = sample.control, destination.folder = file.path(cna.path), Custom.bed=bed, reference.folder = file.path(cna.path, "hg19_400kb"), bp.param = bp.param, bpWidth=1000)

plotCNA(destination.folder = file.path(cna.path))

#load segmented cna profile data:
load(paste0(cna.path,"/","CNAprofiles/segment.Rdata")) # loads a variable called segment.CNA.object

Chromosome=segment.CNA.object$output$chrom
Start=floor(segment.CNA.object$output$loc.start)
End=floor(segment.CNA.object$output$loc.end)
Feature=paste0(Chromosome,":",Start,"-",End)
Mean=segment.CNA.object$output$seg.mean

seg=data.frame(Chromosome=Chromosome, Start=Start, End=End, Feature=Feature, Mean=Mean)
colnames(seg)[5]=name #should make IGV display the sample name properly

header="#track viewLimits=-3:3 graphType=heatmap color=255,0,0"
file.name=paste0(file.path(cna.path),"/",name,".igv")

sink(file.name)
cat(header,"\n")
sink()
write.table(seg, file.name, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, append=TRUE)
