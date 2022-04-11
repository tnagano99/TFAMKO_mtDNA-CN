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
setwd(paste(baseDir, "/results/plots/pathview", sep = ""))

############ KEGG pathview #####################
# KEGG pathway identifiers
# 04080 Neuro 04727 GABA synapse 05033 nicotine 05032 morphine 04974 protein digest absorb 04723 endocannabinoid

# differentially methylated probes
dmps <- read.csv(paste(baseDir, "/results/data/DMPs_anno.csv", sep=""))
dmps <- filter(dmps, pval < 1e-7)

ann <- .getFlatAnnotation("EPIC")
ann$X <- rownames(ann)

df <- merge(dmps, ann, by="X")
df <- df %>% arrange(pval)

geneData <- as.numeric(df$beta)
names(geneData) <- df$entrezid
pathways <- pathview(gene.data = geneData, pathway.id = "04723", species = "hsa", gene.idtype="entrez", out.suffix = "DMP") 

# differentially expressed genes
# 04080 Neuro 04727 GABA synapse 05033 nicotine 05032 morphine 04974 protein digest absorb 04723 endocannabinoid
df <- read.csv(paste(baseDir, "/results/data/EdgeR_RNA_sig_genes.csv", sep=""), header=T) 
names(df)[names(df) == "X"] <- "ensembl_gene_id"

geneData <- as.numeric(df$logFC)
names(geneData) <- df$ensembl_gene_id

pathways <- pathview(gene.data = geneData, pathway.id = "04723", species = "hsa", gene.idtype="ENSEMBL", out.suffix = "RNA")

# differentially methylated region
df1 <- read.csv(paste(baseDir, "/results/data/DMRS_anno.csv", sep=""), header=T) # DMRcate results
df1 <- filter(df1, Fisher < 1.17e-5)
df1 <- separate_rows(df1, overlapping.genes, sep=", ", convert = TRUE)

geneData <- as.numeric(df1$meandiff * 10)
names(geneData) <- df1$overlapping.genes

pathways <- pathview(gene.data = geneData, pathway.id = "04974", species = "hsa", gene.idtype="SYMBOL", out.suffix = "DMR")

ensembl <- mapIds(org.Hs.eg.db, as.character(df1$overlapping.genes), "ENSEMBL", "SYMBOL")
ensembl <- as.vector(unlist(ensembl))
pathways <- pathview(gene.data = ensembl, pathway.id = "04727", species = "hsa", gene.idtype="ENSEMBL", out.suffix = "DMR1")

#Extract genes for KEGG
# z <- getMappedEntrezIDs(sigCpGs, allCpGs, array.type="EPIC")
# keggs <- select(org.Hs.eg.db, z[[1]], "PATH")
# firstkeggs <- subset(keggs, PATH %in% "04080")
# firstkeggs$SYMBOL <- mapIds(org.Hs.eg.db, as.character(firstkeggs$ENTREZID), "SYMBOL","ENTREZID")