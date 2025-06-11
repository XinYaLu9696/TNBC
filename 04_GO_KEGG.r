dev.off()
rm(list=ls())
dir()

Sys.setenv(LANGUAGE = "en") options(stringsAsFactors = FALSE) 
options(clusterProfiler.download.method = "wininet")

library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggThemeAssist)
library(stringr)

gene1<-read.table("../data/up.txt",header=F)
gene2<-read.table("../data/down.txt",header=F)
symbol1=as.character(gene1[,1])
symbol2=as.character(gene2[,1])
eg1 = bitr(symbol1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg2 = bitr(symbol2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
cp = list(Up=eg1$ENTREZID, Down=eg2$ENTREZID)
xy<-compareCluster(cp,fun = "enrichGO",OrgDb="org.Hs.eg.db",keyType="ENTREZID",ont="ALL",
                   pAdjustMethod = "none",pvalueCutoff = 0.05,qvalueCutoff=1,readable=FALSE,pool=TRUE)
aa<-dotplot(xy,split="ONTOLOGY",showCategory=5)+theme_bw()+ 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  scale_color_continuous(low="darkorchid3",high = "sienna2")+facet_grid(ONTOLOGY~., scale='free')
aa+scale_y_discrete(labels=function(y) str_wrap(y,width=60))#use to cut the description to two line
ekk_GO <- setReadable(xy,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
enrich_GO <- as.data.frame(ekk_GO)
head(enrich_GO)
write.csv(enrich_GO,"enrich_GO.csv",quote=F,row.names = F)

xx <- compareCluster(cp, fun="enrichKEGG", organism="hsa", pvalueCutoff=0.05,pAdjustMethod = "none",
                     qvalueCutoff = 1)
bb<-dotplot(xx,showCategory=20)+theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  scale_color_continuous(low="darkorchid3",high = "sienna2")
bb+scale_y_discrete(labels=function(y) str_wrap(y,width=80))
ekk_KEGG <- setReadable(xx,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
enrich_KEGG <- as.data.frame(ekk_KEGG)
head(enrich_KEGG)
write.csv(enrich_KEGG,"enrich_KEGG.csv",quote=F,row.names = F)

