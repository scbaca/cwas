# reformat covariate file to format plink is expecting
args <- commandArgs( trailingOnly = TRUE )
mat = args[1]
out = args[2]

covar=read.table(mat,head=T)
covar=as.data.frame(t(covar[,2:ncol(covar)]))
covar$fam="0"
covar$sample=rownames(covar)
covar = covar[,c((ncol(covar)-1):ncol(covar),1:(ncol(covar)-2))]
write.table(covar, out, quote=F,sep=" ",row.names=F)

