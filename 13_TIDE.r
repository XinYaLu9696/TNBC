dev.off()
rm(list=ls())
dir()
gc()

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
options(clusterProfiler.download.method = "wininet")

library(reshape2)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(scales)
library(limma)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(data.table)
library(stringr)
library(pheatmap)
library(ggsignif)
library(patchwork)

pdata = read.csv("../data/TCGA_FDCSP_TNBC_meta_ok.csv",sep = ",", 
                 row.names = 1, 
                 header = T,check.names = F)
ann = as.data.frame(pdata$group)
rownames(ann) = rownames(pdata)
colnames(ann) = "ImmClust"
head(ann)

TIDE.res <- read.csv("/data/TIDE_TNBC_FDCSP.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
ann$TIDE <- TIDE.res[rownames(ann),"Responder"]
print(table(ann$TIDE,ann$ImmClust))
TIDE.res$cluster = ann[rownames(TIDE.res),"ImmClust"]
res = TIDE.res
NR <- c(38,49)
R <- c(19,7)
dat <- data.frame(NR,R)
rownames(dat) <- c("low","high")
chisq.test(dat)

res = arrange(res,desc(TIDE))
p1 = ggplot(res,aes(x = 1:nrow(res),
                    y = TIDE,
                    fill = Responder))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#e04030","#6cb8d2"))+
  xlab("patient")+
  ylab("TIDE value")+
  #annotate("text", x = 40, y = -1, label = label,size = 10)+
  theme_bw()+
  theme(legend.position = "none") 
p1

res2 = unlist(res)
res2 = as.data.frame(res2)
sapply(res[c("cluster", "Responder")], class)
result <- dplyr::count(res, cluster, Responder)
dat=dplyr::count(res,cluster,Responder)
dat=dat%>%group_by(cluster)%>%
  summarise(Responder=Responder,n=n/sum(n))
dat$Responder=factor(dat$Responder,levels=c("False","True"))
dat

p2=ggplot(data=dat)+
  geom_bar(aes(x=cluster,y=n,
               fill=Responder),
           stat="identity")+
  scale_fill_manual(values=c("#e04030","#6cb8d2"))+
  geom_label(aes(x=cluster,y=n,
                 label=scales::percent(n),
                 fill=Responder),
             color="white",
             size=4,label.size=0,
             show.legend = FALSE,
             position=position_fill(vjust=0.5))+
  ylab("Percentage")+
  theme_minimal()+
  guides(fill = guide_legend(title = "Responder")) 

p1+p2+plot_layout(widths=c(3,2),guides="collect")

res$cluster <- factor(res$cluster,levels = c("low","high"))
png("TIDE_web.png",width = 1200,height = 1200,res = 300)
ggplot(data=res,aes(x=cluster,y=TIDE,colour = cluster))+ 
  geom_violin(
    alpha = 0.8, 
    scale = 'width',
    trim = TRUE)+ 
  geom_boxplot(mapping=aes(x=cluster,y=TIDE,colour=cluster,fill=cluster), 
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=cluster,y=TIDE,colour=cluster), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("low","high"), 
                    values =c("#e04030","#6cb8d2"))+
  scale_color_manual(limits=c("low","high"), 
                     values=c("#e04030","#6cb8d2"))+ 
  geom_signif(mapping=aes(x=cluster,y=TIDE),
              comparisons = list(c("low","high")), 
              map_signif_level=T, 
              tip_length=c(0,0),
              y_position = c(1.5), 
              size=1, 
              textsize = 4, 
              test = "wilcox.test", 
              color = "black")+ 
  theme_bw()+
  guides(fill = guide_legend(title = "cluster"), 
         color = guide_legend(title = "cluster"))+ 
  labs(title = "",  
       x="",y= "TIDE value") 
