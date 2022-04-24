library(limma)
library(biomaRt)
library(edgeR)
library(data.table)
library(dplyr)

############################################################################################

baseDir <- "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN"
setwd(baseDir)

# EdgeR Results
df <- read.csv("./results/data/EdgeR_RNA_all_genes.csv", header=T) # Likelihood test results
names(df)[names(df) == "X"] <- "ensembl_gene_id"

# use biomaRt to find mapping info between ensembl and entrez gene IDS
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) # if unresponsive run: httr::set_config(httr::config(ssl_verifypeer = FALSE))
genes <- getBM(
  attributes=c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
  mart=mart,
  useCache = FALSE)

# drop 3 genes due to missing values in biomart
# EdgeR Filtering and Joing with entrez IDs
df <- merge(df, genes, by = "ensembl_gene_id")
df <- df[!is.na(df$entrezgene_id),]
all_genes <- df$entrezgene_id
df_sig <- filter(df, PValue < 3.53e-6)
df_sig <- df_sig[(df_sig$logFC > 2) | (df_sig$logFC < -2),]
df_sig <- df_sig[!is.na(df_sig$entrezgene_id),]
sig_genes <- df_sig$entrezgene_id

# GO Enrichment
go_rna <- goana(sig_genes, universe = all_genes)
go_rna <- go_rna %>% arrange(P.DE)
go_rna$GO_label <- rownames(go_rna)
write.csv(go_rna, "./results/data/GO_RNA_LRT_EDGER.csv")

# KEGG Enrichment
kegg_rna <- kegga(sig_genes, universe = all_genes)
kegg_rna <- kegg_rna %>% arrange(P.DE)
kegg_rna$KEGG_label <- rownames(kegg_rna)
write.csv(kegg_rna, "./results/data/KEGG_RNA_LRT_EDGER.csv")


################### Old sleuth data ###################
# Sleuth Results Older results no longer used
# df <- read.csv("SleuthWaldTestResults.csv", header=T) # wald test results
# df <- read.csv("SleuthAllGenesAnnotatedRNASeqResultsGeneWise_cleaned.csv", header=T) # Likelihood test results
# names(df)[names(df) == "target_id"] <- "ensembl_gene_id"

# Sleuth Filtering
# merge data with entrez gene ids
# df_anno <- merge(df, genes, by = "ensembl_gene_id")
# df_sig <- filter(df_anno, pval < 3.59e-6) # 3.590149e-06