dev.off()
rm(list=ls())
dir()
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
library(ggplot2)
library(plyr)
library(ggord)
library(yyplot)
library(org.Hs.eg.db)

expr_df = read.csv('ensembl_matrix_TPM.csv', header = T,row.names=1, sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE)  
meta_df = read.csv('meta_data.csv',row.names=1, header = T,sep = ',',quote = '', check.names = FALSE)
expr_df[1:3,1:4]
expr_df = expr_df[rowMeans(expr_df)>10,] 
symbol = as.character(rownames(expr_df))
eg1 = bitr(symbol, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
expr_df$ENSEMBL = rownames(expr_df)
Expr = merge(eg1, expr_df, by="ENSEMBL")
Expr$ENSEMBL = NULL
Expr = avereps(Expr[,-1],ID = Expr$SYMBOL)
Expr = as.data.frame(Expr)
expr_df = as.data.frame(t(Expr))
head(meta_df, n=3)
pca.results <- prcomp(expr_df, center = TRUE, scale. = TRUE)
mycol <- c("#223D6C","#D20A13","#088247","#FFD121","#11AA4D","#58CDD9","#7A142C",
           "#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
ggord(pca.results, grp_in = meta_df$group, repel=TRUE,
      xlims = c(-110,110),
      ylims = c(-60,110),
      ellipse = FALSE, 
      size = 2, 
      alpha=0.5, 
      cols = mycol[1:length(unique(meta_df$group))],
      arrow = NULL,txt = NULL) +
  theme(panel.grid =element_blank()) 


