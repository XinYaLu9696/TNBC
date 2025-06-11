dev.off()
rm(list=ls())
dir()
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)


data = read.csv('../data/Expr_fpkm_TCGA_ok.csv', 
                header = T,row.names=1, sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE)  
data = as.data.frame(t(data))
pdata = read.csv('../data/meta_merge_OK.csv', 
                 header = T,row.names=1, sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE)  
data = data[rownames(data)%in%rownames(pdata),]
pdata = pdata[rownames(pdata)%in%rownames(data),]
data = data[rownames(pdata),]
all(rownames(data)==rownames(pdata))
data_need = data[,c("FDCSP","FOXC1","LCN2","MMP7","KRT16","KRT23","KRT6A","KRT6B")]
plot_data = cbind(pdata,data_need)
plot_data = plot_data[,c(2,3)]
table(plot_data$group)
colnames(plot_data) = c("group","Retive_Abundance")
plot_data$Retive_Abundance = log2(plot_data$Retive_Abundance+1)
plot_data$group <- factor(plot_data$group,levels = c("Normal","no_TNBC","TNBC"))
p<- ggplot(data=plot_data)+ 
  geom_boxplot(mapping=aes(x=group,y=(Retive_Abundance),colour = group ), 
               alpha = 0.5,
               size=0.75,
               width = 0.6)+ 
  geom_jitter(mapping=aes(x=group,y=Retive_Abundance,colour = group), 
              width = 0.2, alpha = 1,size=1)+
  scale_x_discrete(limits = c("Normal","no_TNBC","TNBC"))+
  scale_color_manual(limits=c("Normal","no_TNBC","TNBC"), 
                     values=c("#45B39D", "#F5B041","#A93226","darkred","#922927","#922927"))+ 
  geom_signif(mapping=aes(x=group,y=Retive_Abundance), 
              comparisons = list(c("Normal", "no_TNBC"), 
                                  c("no_TNBC", "TNBC"),
                                  c("Normal", "TNBC") 
                                 ),
              map_signif_level=T, 
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(17,18,19), 
              size=0.75, 
              textsize = 4, 
              test = wilcox.test)+ 
  theme_classic( 
    base_line_size = 0.75 
  )+
  labs(title="FDCSP",x="",y="log2(FPKM+1)")+ 
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "plain", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color="black", 
                                    size=15, 
                                    face="plain"),
        legend.text = element_text(color="black", 
                                   size = 10, 
                                   face = "plain"),
        axis.text.x = element_text(size = 13,  
                                   color = "black", 
                                   face = "plain", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0), 
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "plain", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )
p
