# compare average Z-squared of cwas associations across different gwas studies
library(ggplot2)
library(dplyr)

args=commandArgs(trailingOnly = TRUE)

#tmp
#args=c('analysis/fusion/ProstateCancer_Meta_Schumacher2018.nodup/merged/cv/cwas.cv.best.txt','analysis/fusion/UKB_460K.biochemistry_Testosterone_Male/merged/cv/cwas.cv.best.txt','analysis/fusion/UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED/merged/cv/cwas.cv.best.txt','analysis/fusion/UKB_460K.disease_PSORIASIS/merged/cv/cwas.cv.best.txt','analysis/fusion/UKB_460K.disease_DERMATOLOGY/merged/cv/cwas.cv.best.txt','analysis/fusion/UKB_460K.disease_AID_ALL/merged/cv/cwas.cv.best.txt','analysis/fusion/UKB_460K.disease_ASTHMA_DIAGNOSED/merged/cv/cwas.cv.best.txt','analysis/fusion/UKB_460K.disease_CARDIOVASCULAR/merged/cv/cwas.cv.best.txt','analysis/fusion/UKB_460K.disease_HI_CHOL_SELF_REP/merged/cv/cwas.cv.best.txt','analysis/fusion/UKB_460K.disease_HYPERTENSION_DIAGNOSED/merged/cv/cwas.cv.best.txt','analysis/fusion/UKB_460K.disease_RESPIRATORY_ENT/merged/cv/cwas.cv.best.txt','analysis/fusion/UKB_460K.disease_T2D/merged/cv/cwas.cv.best.txt','analysis/fusion/UKB_460K.disease_HYPOTHYROIDISM_SELF_REP/merged/cv/cwas.cv.best.txt')

plotfile="analysis/fusion/compare_gwas/compare.gwas.pdf"
n=length(args)
names=gsub("analysis/fusion/","",args)
names=gsub("/merged/cv/cwas.cv.best.txt","",names)

zs=lapply(seq_along(args), function(x){df=read.table(args[x], header=T); df$gwas=names[x]; df}) %>% 
  do.call("rbind", .)

zs$abs_z=abs(zs$TWAS.Z)

#limit to top 100
zs = zs %>% group_by(gwas) %>% 
  top_n(wt = abs_z, n = 100) %>% 
  mutate(mean=mean(abs_z)) %>% 
  as.data.frame()

names=as.factor(names)
ggplot(zs, aes(y=abs_z,x=reorder(gwas, -mean))) +
  geom_jitter(width = 0.25, col="darkgray") + 
  geom_boxplot(outlier.shape=NA) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  xlab("") + 
  ylab("|z-score|") +

ggsave(plotfile, height=5, width=4)
