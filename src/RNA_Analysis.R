library(limma)
library(biomaRt)
library(edgeR)
library(data.table)
library(dplyr)

# setwd("../results/data")

# df <- read.csv("SleuthGeneNormalized_TPM.csv", header=T)
# row.names(df) <- df$X
# df$X <- NULL
# colnames(df) <- c("TFAM_KO_Clone_11_1", "TFAM_KO_Clone_11_2", "TFAM_KO_Clone_5_1", "TFAM_KO_Clone_5_2", "TFAM_KO_Clone_6_1", "TFAM_KO_Clone_6_2", "HEK293T_NC.1_1", "HEK293T_NC.1_2", "HEK293T_NC.2_1", "HEK293T_NC.2_2", "HEK293T_NC.3_1", "HEK293T_NC.3_2")
# df <- df[, colnames(EPIC)]

# #Remove genes with TPM < 0.5 in >49% of samples
# df$TPM05 <- rowSums(df<0.5)
# df <- subset(df, df$TPM05 < 7)
# df$TPM05 <- NULL

# # edgeR
# # specify conditions and design matrix for limma
# type <- factor(c("Experiment", "Control", "Experiment", "Control", "Experiment", "Control", "Control", "Control", "Control", "Experiment", "Experiment", "Experiment"))
# batch <- factor(c(rep(203751390048, 6), rep(201114530008, 6)))


############################################################################################

setwd("results/data")
# df <- read.csv("SleuthWaldTestResults.csv", header=T) # wald test results
df <- read.csv("SleuthAllGenesAnnotatedRNASeqResultsGeneWise_cleaned.csv", header=T) # Likelihood test results
names(df)[names(df) == "target_id"] <- "ensembl_gene_id"

# use biomaRt to find mapping info between ensembl and entrez gene IDS
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) # if unresponsive run: httr::set_config(httr::config(ssl_verifypeer = FALSE))
genes <- getBM(
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  mart=mart,
  useCache = FALSE)

# merge data with entrez gene ids
df_anno <- merge(df, genes, by = "ensembl_gene_id")
df_sig <- filter(df_anno, pval < 3.59e-6) # 3.590149e-06
# df_sig <- filter(df_anno, qval < 0.05)
df_sig <- df_sig[!is.na(df_sig$entrezgene_id),]
sig_genes <- df_sig$entrezgene_id

go_rna <- goana(sig_genes, universe = genes$entrezgene_id)
go_rna <- go_rna %>% arrange(P.DE)
go_rna$GO_label <- rownames(go_rna)
write.csv(go_rna, "GO_RNA_LRT.csv")

kegg_rna <- kegga(sig_genes, universe = genes$entrezgene_id)
kegg_rna <- kegg_rna %>% arrange(P.DE)
kegg_rna$KEGG_label <- rownames(kegg_rna)
write.csv(kegg_rna, "KEGG_RNA_LRT.csv")