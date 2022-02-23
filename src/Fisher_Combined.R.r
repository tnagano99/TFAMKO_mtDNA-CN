library(edgeR)
library(limma)
library(biomaRt)
library(metap)
library(dplyr)
library(data.table)

# set working directory to be folder containing root
baseDir <- "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN"
setwd(baseDir)

# read in csv data for GO and KEGG
go <- as.data.frame(read.csv(paste(baseDir, "/results/data/GO_Cont.csv", sep = "")))
# go <- as.data.frame(read.csv(paste(baseDir, "/results/data/GO_DMRcate_pval.csv", sep = "")))
go_rna <- as.data.frame(read.csv(paste(baseDir, "/results/data/GO_RNA_LRT_new.csv", sep = "")))

kegg <- as.data.frame(read.csv(paste(baseDir, "/results/data/KEGG_Cont.csv", sep = "")))
# kegg <- as.data.frame(read.csv(paste(baseDir, "/results/data/KEGG_DMRcate_pval.csv", sep = "")))
kegg_rna <- as.data.frame(read.csv(paste(baseDir, "/results/data/KEGG_RNA_LRT_new.csv", sep = "")))

# Combine GO and KEGG results using fisher method
# rename columns to combine p values and merge data
colnames(go)[colnames(go) == "P.DE"] <- "P.DE_Meth"
colnames(go)[colnames(go) == "X"] <- "GO_label"
colnames(go_rna)[colnames(go_rna) == "P.DE"] <- "P.DE_RNA"
colnames(kegg)[colnames(kegg) == "P.DE"] <- "P.DE_Meth"
colnames(kegg)[colnames(kegg) == "X"] <- "KEGG_label"
colnames(kegg_rna)[colnames(kegg_rna) == "P.DE"] <- "P.DE_RNA"

go_rna_sub <- go_rna[c("P.DE_RNA", "GO_label")]
kegg_rna_sub <- kegg_rna[c("P.DE_RNA", "KEGG_label")]

go_all <- merge(go, go_rna_sub, by= "GO_label")
kegg_all <- merge(kegg, kegg_rna_sub, by= "KEGG_label")

# combine p-value using fisher method
go_comb <- numeric(dim(go_all)[1])
for (i in 1:dim(go_all)[1]){
  p_vals = c(go_all$P.DE_Meth[i], go_all$P.DE_RNA[i])
  fisher <- sumlog(p_vals)
  go_comb[i] <- fisher$p
}

kegg_comb <- numeric(dim(kegg_all)[1])
for (i in 1:dim(kegg_all)[1]){
  p_vals = c(kegg_all$P.DE_Meth[i], kegg_all$P.DE_RNA[i])
  fisher <- sumlog(p_vals)
  kegg_comb[i] <- fisher$p
}

go_all$Fisher <- go_comb
go_all <- go_all %>% arrange(Fisher)
write.csv(go_all, paste(baseDir, "/results/data/GO_Probe_RNA_Fisher.csv", sep=""))

kegg_all$Fisher <- kegg_comb
kegg_all <- kegg_all %>% arrange(Fisher)
write.csv(kegg_all, paste(baseDir, "/results/data/KEGG_Probe_RNA_Fisher.csv", sep=""))