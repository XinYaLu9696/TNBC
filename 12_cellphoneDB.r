cellphonedb method statistical_analysis meta_Normal.txt count_Normal.txt --database cellphone.db --counts-data=gene_name --output-path ./out
cellphonedb method statistical_analysis meta_ER.txt count_ER.txt --database cellphone.db --counts-data=gene_name --output-path ./out
cellphonedb method statistical_analysis meta_HER2.txt count_HER2.txt --database cellphone.db --counts-data=gene_name --output-path ./out
cellphonedb method statistical_analysis meta_TNBC.txt count_TNBC.txt --database cellphone.db --counts-data=gene_name --output-path ./out

#====================================================================================================================================
dev.off()
rm(list=ls())
gc()
dir()
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
library(ktplots)
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(data.table)
library(dplyr)
library(readr)
mypvals <- read.table("pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
mymeans <- read.table("means.txt",header = T,sep = "\t",stringsAsFactors = F) 
kp = grepl(pattern = "THBS2", colnames(mypvals)) & grepl(pattern = "capillary", colnames(mypvals)) 
table(kp)
pos = (1:ncol(mypvals))[kp] 
choose_pvalues <- mypvals[,c(c(1,5,6,8,9),pos)]
choose_means <- mymeans[,c(c(1,5,6,8,9),pos)]
logi <- apply(choose_pvalues[,6:ncol(choose_pvalues)]<0.05, 1, sum) 
choose_pvalues <- choose_pvalues[logi>=1,]
logi1 <- choose_pvalues$gene_a != ""
logi2 <- choose_pvalues$gene_b != ""
logi <- logi1 & logi2
choose_pvalues <- choose_pvalues[logi,]
choose_means <- choose_means[choose_means$id_cp_interaction %in% choose_pvalues$id_cp_interaction,]


meansdf <- choose_means %>% reshape2::melt()
meansdf <- data.frame(interacting_pair = paste0(meansdf$gene_a,"_",meansdf$gene_b),
                      CC = meansdf$variable,
                      means = meansdf$value)
pvalsdf <- choose_pvalues %>% reshape2::melt()
pvalsdf <- data.frame(interacting_pair = paste0(pvalsdf$gene_a,"_",pvalsdf$gene_b),
                      CC = pvalsdf$variable,
                      pvals = pvalsdf$value)
pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")

summary((filter(pldf,means >0))$means)
head(pldf)
pcc =  pldf%>% filter(means >0) %>% 
  ggplot(aes(CC.x,interacting_pair.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals+0.0001)) ) +
  #scale_size_continuous(range = c(1,3))+
  scale_color_gradient2(high="red",mid = "grey",low ="darkblue",midpoint = 0.9  )+ 
  #scale_color_gradient2(high="Firebrick",mid = "Khaki1",low ="SlateBlue3",midpoint = 0.6  )+ 
  theme_bw()+ 
  # scale_color_manual(values = rainbow(100))+
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 0))+
  labs(x="Cell interacting",y="interacting pair")
pcc

