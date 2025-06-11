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
library(DOSE)
library(stringr)

phenotype_file = read.table('Clinical BCR XML.merge.txt',header = T,sep = '\t',quote = '')
table(phenotype_file$breast_carcinoma_estrogen_receptor_status) 
table(phenotype_file$breast_carcinoma_progesterone_receptor_status) 
table(phenotype_file$lab_proc_her2_neu_immunohistochemistry_receptor_status) 
colnames_num <- grep('receptor_status',colnames(phenotype_file))
colnames_num
phenotype_colnames <- colnames(phenotype_file)[colnames_num]
phenotype_colnames
eph <- phenotype_file[,colnames_num[c(1,2,12)]]
write.csv(eph, file = "phenotype_file_OK.csv", quote=F, row.names=T)
tnbc_rownum <- apply(eph, 1, function(x) sum(x =='Negative'))
tnbc_sample <- phenotype_file[tnbc_rownum == 3, 1]
tnbc_sample
no_tnbc_sample <- phenotype_file[tnbc_rownum != 3, 1]
write.csv(tnbc_sample, file = "tnbc_sample.csv", quote=F, row.names=F)
write.csv(no_tnbc_sample, file = "no_tnbc_sample.csv", quote=F, row.names=F)
