dev.off()
rm(list=ls())
dir()
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

library(limma)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(data.table)
library(topGO)
library(GSEABase)
library(stringr)

a = read.csv('count_matrix.csv', header = T, sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE) 
a = avereps(a[,-1],ID = a$GENE) 
a = as.data.frame(a)
Expr = a
pdata = read.csv('meta_data.csv',row.names=1,header = T,sep = ',',quote = '',check.names = FALSE)
pdata <- as.data.frame(pdata)
pdata[1:3,]

pdata = pdata[colnames(Expr),]
colnames(Expr) == rownames(pdata)

aa = 0
ak <- function(x) {
  if (as.numeric(x) > aa) {
    x = as.factor("high")
    return(x)
  } else {
    x = as.factor("low")
    return(x)
  }
}

pp <- data.frame(sapply(pdata$TNBC, function(x) ak(x)))
rownames(pp) = rownames(pdata)
colnames(pp) = "type" 

#Step1
#limma
group_list = as.character(pp[, 1])
table(group_list)
# Multidimensional scaling (MDS) plot
plotMDS(Expr, col = as.numeric(group_list))

# Step2
design <- model.matrix(~0+factor(group_list))
head(design)
levels(factor(group_list))
colnames(design) = levels(factor(group_list))
rownames(design) = rownames(pdata)

# Step3
DGElist = DGEList(counts = Expr, group = group_list)
keep_gene = rowSums( cpm(DGElist) > 1 ) >= 10
table(keep_gene)
DGElist = DGElist[keep_gene, , keep.lib.sizes = FALSE]

# Step4
DGElist = calcNormFactors(DGElist, method = 'TMM') 
v = voom(DGElist, design, plot = TRUE, normalize = "quantile")   
fit = lmFit(v, design)
cont.matrix = makeContrasts(contrasts = c('high-low'), levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)
plotSA(fit2, main="Final model: Mean-variance trend")

# Step5
nrDEG_limma_voom = topTable(fit2, coef = 'high-low', n = Inf, adjust="BH")
nrDEG_limma_voom = na.omit(nrDEG_limma_voom)
head(nrDEG_limma_voom)

write.csv(nrDEG_limma_voom, file = "DEG_matrix_TNBCvsother.csv", quote=F, row.names=T)
