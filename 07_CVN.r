dev.off()
rm(list=ls())
dir()
gc()
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

library(maftools)
library(TCGAmutations)
library(limma)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(data.table)
library(stringr)

tcga_mc3 = tcga_load(study = "BRCA") 
save(tcga_mc3, file = "TCGA_BRCA_MC3_maf.Rdata")
load("TCGA_BRCA_MC3_maf.Rdata")
d1 <- as.data.frame(tcga_mc3@data)
write.csv(d1, file="tcga_brca_maf.csv", quote=F, row.names=T)  
meta = read.csv('../data/TCGA_FDCSP_TNBC_meta_ok.csv', header = T,row.names=1, sep = ',', quote = '', 
                    stringsAsFactors = FALSE, check.names = FALSE)  
meta$sample = substr(rownames(meta),1,12)
meta_TNBC = meta[meta$group == "low",]
d1_TNBC = d1[d1$Tumor_Sample_Barcode_min %in% meta_TNBC$sample,]
ab_high = d1_TNBC[,c("Tumor_Sample_Barcode","Matched_Norm_Sample_Barcode")]
ab_high = unique(ab_high)
meta_TNBC = meta_TNBC[ab_high$Tumor_Sample_Bar,]
pdata_high_merge = cbind(meta_TNBC,ab_high)
meta_no = meta[meta$group == "high",]
d1_no = d1[d1$Tumor_Sample_Barcode_min %in% meta_no$sample,]
ab_low = d1_no[,c("Tumor_Sample_Barcode","Matched_Norm_Sample_Barcode")]
ab_low = unique(ab_low)
meta_no = meta_no[ab_low$Tumor_Sample_Bar,]
pdata_low_merge = cbind(meta_no,ab_low)
laml = read.maf(maf = d1_TNBC, clinicalData = pdata_high_merge) 
#laml = read.maf(maf = d1_no, clinicalData = pdata_low_merge) 
getSampleSummary(laml)
#Shows gene summary.
getGeneSummary(laml)
#Shows all fields in MAF
getFields(laml)
#shows clinical data associated with samples
getClinicalData(laml)
#Writes maf summary to an output file with basename laml.
#write.mafSummary(maf = laml, basename = 'laml')
#plotmafSummary
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', 
               dashboard = TRUE, titvRaw = FALSE)
laml = tcga_mc3
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
#oncoplot for top ten mutated genes.
oncoplot(maf = laml, top = 20)


