dev.off()
rm(list=ls())
dir()
gc()

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
options(clusterProfiler.download.method = "wininet")

set.seed(1)
library(tidyverse)
library(pheatmap)
library(CIBERSORT)
library(limma)

a = read.csv('../data/Expr_fpkm_TCGA_ok.csv', header = T,sep = ',', quote = '', 
             stringsAsFactors = FALSE, check.names = FALSE)  
exp  = avereps(a[,-1],ID = a[,1]) 
exp  = as.data.frame(exp)
exp2 = exp[rowSums(exp>=1)>=10,]
Expr_count = rownames_to_column(exp2 )

f = "ciber_mouse_process.Rdata"
if(!file.exists(f)){
  #devtools:: install_github ("Moonerss/CIBERSORT")
  lm22f = system.file("extdata", "LM22.txt", package = "CIBERSORT")
  TME.results = cibersort(lm22f, 
                          "exp2.txt" , 
                          perm = 1000, 
                          QN = T)
  save(TME.results,file = f)
}


load(f)
TME.results[1:4,1:4]
re <- TME.results[,-(23:25)]
a = read.csv('../data/TCGA_FDCSP_TNBC_meta_ok.csv', header = T,row.names=1, 
             sep = ',', quote = '', 
             stringsAsFactors = FALSE, check.names = FALSE)  
re = re[rownames(a),]
library(pheatmap)
k <- apply(re,2,function(x) {sum(x == 0) < nrow(TME.results)/2})
table(k)
re2 <- as.data.frame(t(re[,k]))
all(colnames(re2) == rownames(a))
Group = a$group
an = data.frame(group = Group,
                row.names = colnames(re2))

library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
dat <- re %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  mutate(group = Group) %>%
  gather(key = Cell_type,value = Proportion,-Sample,-group) %>%
  arrange(group)

dat$Sample = factor(dat$Sample,ordered = T,levels = unique(dat$Sample)) 
dat2 = data.frame(a = 1:ncol(re2),
                  b = 1,
                  group = sort(Group))

p1 = ggplot(dat2,aes(x = a, y = b)) +
  geom_tile(aes(fill = group)) +
  scale_fill_manual(values = mypalette(22)[1:length(unique(Group))]) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(fill = "Group")
p1
p2 = ggplot(dat,aes(Sample, Proportion,fill = Cell_type)) +
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))
p2
library(patchwork)
p1 / p2 + plot_layout(heights = c(1,10),guides = "collect" ) &
  theme(legend.position = "bottom")

library(vioplot)         
normal= 57
tumor=  56
vioplotdata = as.data.frame((re))
all(rownames(vioplotdata) == rownames(Group))
vioplotdata = vioplotdata[,which(colSums(vioplotdata) > 0)] 
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(vioplotdata))
y=c(1:ncol(vioplotdata))
plot(x,y,
     xlim=c(0,(ncol(vioplotdata)-1)*3),
     #xlim=c(0,57)
     ylim=c(min(vioplotdata),max(vioplotdata)+0.02),
     legend=TRUE,
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")
for(i in 1:ncol(vioplotdata)){
  normalData=vioplotdata[1:normal,i]
  tumorData=vioplotdata[(normal+1):(normal+tumor),i]
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = ("#5390d9"))
  #vioplot(normalData,at=3*(i-1),lty=1,add = T,
   #       col = ("#a3de83"))
  #vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col =  ("#dc2f02"))
  vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,
          col =  ("#a3de83"))
  wilcoxTest=wilcox.test(normalData,tumorData)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData))
  #lines(c(x=3*(i-1)+0.3,ax=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5,y=mx+0.02,labels=ifelse(p<0.05,paste0("p<0.05"),paste0("p=",p)),cex = 0.8)
  text(seq(1,(ncol(vioplotdata)-1)*3+1,3),-0.05,xpd = NA,labels=colnames(vioplotdata),cex = 1,srt = 45,pos=2)
}
