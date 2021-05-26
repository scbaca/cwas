#plot stats for enrichment of eQTLs in specified peak sets. NOTE: qtl files must have "qtl" in the title for this to work.

library(reshape)
library(ggplot2)
library(dplyr)
library(stringr)

args <- commandArgs( trailingOnly = TRUE )
summ_file=args[1]
plot_file=args[2]

d=read.table(summ_file, sep="\t", header=F, colClasses=c("character", rep("NULL",5), "numeric", rep("NULL",2), "numeric", "NULL")) 
colnames(d) = c("dataset", "enrichment", "enrichment.rand")
d$dataset = d$dataset %>% str_replace("qtl/", "") %>% str_replace(".qtl.hg19", "")

d=melt(d)

#order by enrichment:
lev = d$dataset[d$variable=="enrichment.rand"]
lev = lev[order(d$value[d$variable=="enrichment.rand"], decreasing=T)]
d$dataset = factor(d$dataset, levels = lev)

ggplot(d, aes(y=value, x=dataset, fill=variable, width=0.6)) + geom_col(position="dodge") + theme_classic() + ylab("enrichment") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ggtitle("eQTL overlap enrichment")
ggsave(plot_file, height=3, width=4)             
