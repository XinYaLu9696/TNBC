
dev.off()
rm(list=ls())
dir()
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE)
library(pheatmap)
library(RColorBrewer)
expr_df = read.csv('Expr_tpm_ok.csv', header = T,row.names=1, sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE)  
meta_df = read.csv('meta_data.csv',row.names=1,header = T,sep = ',',quote = '', check.names = FALSE)
Expr = as.data.frame(expr_df)
data = read.table("genes.txt")
a = Expr[data$V1,]
bk = unique(c(seq(-1,1, length=100)))
aa = pheatmap((a), 
              scale = "row", 
              breaks = bk,
              clustering_distance_cols = "correlation", 
              clustering_method="average",
              fontsize_row = 10,
              show_colnames = T, 
              show_rownames = T, 
              cluster_cols = F,
              cluster_row = T, 
              color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")[1:7]))(100)
)

