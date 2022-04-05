library(edgeR)
library(limma)
library(metap)
library(dplyr)
library(data.table)

# set working directory to be folder containing root
baseDir <- "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN"
setwd(baseDir)

################## GO/KEGG Fisher Combined ##########################
# swap go and kegg variables for probe and region combined
# read in csv data for Meth and RNA GO results
# go <- as.data.frame(read.csv(paste(baseDir, "/results/data/GO_Cont.csv", sep = ""))) # Probe
go <- as.data.frame(read.csv(paste(baseDir, "/results/data/GO_DMRcate.csv", sep = ""))) # RNA
go_rna <- as.data.frame(read.csv(paste(baseDir, "/results/data/GO_RNA_LRT_EDGER.csv", sep = "")))

# read in csv data for Meth and RNA KEGG results
# kegg <- as.data.frame(read.csv(paste(baseDir, "/results/data/KEGG_Cont.csv", sep = ""))) # Probe
kegg <- as.data.frame(read.csv(paste(baseDir, "/results/data/KEGG_DMRcate.csv", sep = ""))) # RNA
kegg_rna <- as.data.frame(read.csv(paste(baseDir, "/results/data/KEGG_RNA_LRT_EDGER.csv", sep = "")))

# Combine GO and KEGG results using fisher method
# rename columns to combine p values and merge data
colnames(go)[colnames(go) == "P.DE"] <- "P.DE_Meth"
colnames(go)[colnames(go) == "X"] <- "GO_label"
colnames(go_rna)[colnames(go_rna) == "P.DE"] <- "P.DE_RNA"
colnames(kegg)[colnames(kegg) == "P.DE"] <- "P.DE_Meth"
colnames(kegg)[colnames(kegg) == "X"] <- "KEGG_label"
colnames(kegg_rna)[colnames(kegg_rna) == "P.DE"] <- "P.DE_RNA"

# remove entries with p > 0.05 from both meth and RNA
go <- subset(go, P.DE_Meth < 0.05)
go_rna <- subset(go_rna, P.DE_RNA < 0.05)
kegg <- subset(kegg, P.DE_Meth < 0.05)
kegg_rna <- subset(kegg_rna, P.DE_RNA < 0.05)

# subset the RNA results to combine with Meth
go_rna_sub <- go_rna[c("P.DE_RNA", "GO_label")]
kegg_rna_sub <- kegg_rna[c("P.DE_RNA", "KEGG_label")]

# merge RNA and Meth results
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

# add GO/KEGG Fisher combined value to dataframes and save
go_all$Fisher <- go_comb
go_all <- go_all %>% arrange(Fisher)
write.csv(go_all, paste(baseDir, "/results/data/GO_Region_RNA_Fisher_EDGER.csv", sep=""))

kegg_all$Fisher <- kegg_comb
kegg_all <- kegg_all %>% arrange(Fisher)
write.csv(kegg_all, paste(baseDir, "/results/data/KEGG_Region_RNA_Fisher_EDGER.csv", sep=""))

################ TFT/Reactome/MitoCarta Fisher Combined ################

# swap go and kegg variables for probe and region combined
# read in csv data for TFT Meth and RNA results
tft <- as.data.frame(read.csv(paste(baseDir, "/results/data/overrepresentation_tft_meth.csv", sep = ""))) 
tft_rna <- as.data.frame(read.csv(paste(baseDir, "/results/data/overrepresentation_tft_log.csv", sep = "")))

# read in csv data for Reactome Meth and RNA results
# reactome <- as.data.frame(read.csv(paste(baseDir, "/results/data/overrepresentation_reactome_meth.csv", sep = "")))
# reactome_rna <- as.data.frame(read.csv(paste(baseDir, "/results/data/overrepresentation_reactome_log.csv", sep = "")))

# mitocarta no overlapping pathways
# mt <- as.data.frame(read.csv(paste(baseDir, "/results/data/overrepresentation_mt_meth.csv", sep = "")))
# mt_rna <- as.data.frame(read.csv(paste(baseDir, "/results/data/overrepresentation_mt_log.csv", sep = "")))

# Combine GO and KEGG results using fisher method
# rename columns to combine p values and merge data
colnames(tft)[colnames(tft) == "T.test.pval"] <- "P.DE_Meth"
colnames(tft_rna)[colnames(tft_rna) == "T.test.pval"] <- "P.DE_RNA"
# colnames(reactome)[colnames(reactome) == "T.test.pval"] <- "P.DE_Meth"
# colnames(reactome_rna)[colnames(reactome_rna) == "T.test.pval"] <- "P.DE_RNA"
# colnames(mt)[colnames(mt) == "T.test.pval"] <- "P.DE_Meth"
# colnames(mt_rna)[colnames(mt_rna) == "T.test.pval"] <- "P.DE_RNA"

# remove entries with p > 0.05 from both meth and RNA
tft <- subset(tft, P.DE_Meth < 0.05)
tft_rna <- subset(tft_rna, P.DE_RNA < 0.05)
# reactome <- subset(reactome, P.DE_Meth < 0.05)
# reactome_rna <- subset(reactome_rna, P.DE_RNA < 0.05)
# mt <- subset(mt, P.DE_Meth < 0.05)
# mt_rna <- subset(mt_rna, P.DE_RNA < 0.05)

# subset the RNA results to combine with Meth
tft_rna_sub <- tft_rna[c("P.DE_RNA", "Gene.Set.Name")]
# reactome_rna_sub <- reactome_rna[c("P.DE_RNA", "Gene.Set.Name")]
# mt_rna_sub <- mt_rna[c("P.DE_RNA", "Gene.Set.Name")]

# merge RNA and Meth results
tft_all <- merge(tft, tft_rna_sub, by= "Gene.Set.Name")
# reactome_all <- merge(reactome, reactome_rna_sub, by= "Gene.Set.Name")
# mt_all <- merge(mt, mt_rna_sub, by= "Gene.Set.Name")

# combine p-value using fisher method
tft_comb <- numeric(dim(tft_all)[1])
for (i in 1:dim(tft_all)[1]){
  p_vals = c(tft_all$P.DE_Meth[i], tft_all$P.DE_RNA[i])
  fisher <- sumlog(p_vals)
  tft_comb[i] <- fisher$p
}

# reactome_comb <- numeric(dim(reactome_all)[1])
# for (i in 1:dim(reactome_all)[1]){
#   p_vals = c(reactome_all$P.DE_Meth[i], reactome_all$P.DE_RNA[i])
#   fisher <- sumlog(p_vals)
#   reactome_comb[i] <- fisher$p
# }

# mt_comb <- numeric(dim(mt_all)[1])
# for (i in 1:dim(mt_all)[1]){
#   p_vals = c(mt_all$P.DE_Meth[i], mt_all$P.DE_RNA[i])
#   fisher <- sumlog(p_vals)
#   mt_comb[i] <- fisher$p
# }

tft_all$Fisher <- tft_comb
tft_all <- tft_all %>% arrange(Fisher)
write.csv(tft_all, paste(baseDir, "/results/data/TFT_Probe_RNA_Fisher.csv", sep=""))

# reactome_all$Fisher <- reactome_comb
# reactome_all <- reactome_all %>% arrange(Fisher)
# write.csv(reactome_all, paste(baseDir, "/results/data/REACTOME_Probe_RNA_Fisher.csv", sep=""))

# mt_all$Fisher <- mt_comb
# mt_all <- mt_all %>% arrange(Fisher)
# write.csv(mt_all, paste(baseDir, "/results/data/MITOCARTA_Probe_RNA_Fisher.csv", sep=""))