library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(minfi)
library(missMethyl)
library(pathview)
library(biomaRt)
library("org.Hs.eg.db")
# library(qqman)
# library(limma)

# Update baseDir to proper location
baseDir <- "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN"
setwd(baseDir)

############ KEGG pathview #####################
dmps <- read.csv(paste(baseDir, "/results/data/DMPs_anno.csv", sep=""))

# differentially methylated probes
subset <- as.data.frame(dmps[,c("X","GencodeCompV12_Accession")])
expanded <- separate_rows(subset, GencodeCompV12_Accession, sep=";", convert = TRUE)
test <- transpose(as.data.frame(str_split(expanded$GencodeCompV12_Accession, '[.]', n = 2)))
expanded$GencodeCompV12_Accession <- NULL
expanded$ensembl_gene_id <- test[,1]

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) # if unresponsive run: httr::set_config(httr::config(ssl_verifypeer = FALSE))
genes <- getBM(
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  mart=mart,
  useCache = FALSE)

df_anno <- merge(expanded, genes, by = "ensembl_gene_id")

pathways <- pathview(gene.data = expanded$ensembl_gene_id, pathway.id = "05033", species = "hsa", gene.idtype="ENSEMBLTRANS", kegg.dir = "./results/plots", out.suffix = "DMP") # 04080 Neuro 04727 GABA synapse 05033 nicotine

# differentially expressed genes
# df <- read.csv("./results/data/SleuthAllGenesAnnotatedRNASeqResultsGeneWise_cleaned.csv", header=T) # Likelihood test results
df <- read.csv("./results/data/EdgeR_RNA_all_genes.csv", header=T) 
names(df)[names(df) == "X"] <- "ensembl_gene_id"

df <- filter(df, pval < 3.59e-6)
pathways <- pathview(gene.data = df$ensembl_gene_id, pathway.id = "05033", species = "hsa", gene.idtype="ENSEMBL", kegg.dir = "./results/plots", out.suffix = "RNA")

# differentially methylated region
df1 <- read.csv("./results/data/DMRS_anno.csv", header=T) # DMRcate results
df1 <- filter(df1, Fisher < 1.17e-5)
df1 <- separate_rows(df1, overlapping.genes, sep=", ", convert = TRUE)
pathways <- pathview(gene.data = df1$overlapping.genes, pathway.id = "05033", species = "hsa", gene.idtype="SYMBOL", kegg.dir = "./results/plots", out.suffix = "DMR")

df1$ENSEMBL <- mapIds(org.Hs.eg.db, as.character(df1$overlapping.genes), "ENSEMBL", "SYMBOL")
df1 <- as.data.frame(df1)
pathways <- pathview(gene.data = df1$ENSEMBL, pathway.id = "05033", species = "hsa", gene.idtype="ENSEMBL", kegg.dir = "./results/plots", out.suffix = "DMR")

#Extract genes for KEGG
z <- getMappedEntrezIDs(sigCpGs, allCpGs, array.type="EPIC")
keggs <- select(org.Hs.eg.db, z[[1]], "PATH")
firstkeggs <- subset(keggs, PATH %in% "04080")
firstkeggs$SYMBOL <- mapIds(org.Hs.eg.db, as.character(firstkeggs$ENTREZID), "SYMBOL","ENTREZID")