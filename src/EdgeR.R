library(edgeR)
library(dplyr)
library(data.table)
library(biomaRt)

##################### Format Kaillisto output to tracsript-level counts #####################
# getting estimated counts from Kallisto output
baseDir <- "/home/tnagano/projects/def-ccastel/SharedResources/TFAM_KO/FASTQ"

# use the metaData from csv to rename columns
metaData <- read.csv(paste(baseDir, "targetsSleuthRD.csv", sep="/"))
df <- read.csv(paste(baseDir, "kallisto_output1/abundance.tsv", sep="/"), sep = "\t")
df <- df[,c("target_id","est_counts")]
colnames(df)[colnames(df) == "est_counts"] <- metaData$sample[1]

# loop through abundance.tsv files to get estimated_count outputs from kallisto
for (i in 2:12){
    subDir <- paste0("kallisto_output", i)
    tmp <- read.csv(paste(baseDir, subDir,"abundance.tsv", sep="/"), sep = "\t")
    tmp <- tmp[,c("target_id","est_counts")]
    colnames(tmp)[colnames(tmp) == "est_counts"] <- metaData$sample[i]
    df <- merge(df, tmp, by = "target_id")
}

# write.csv(df, "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN/results/data/EdgeRTranscriptEstCounts.csv")
df <- read.csv("/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN/results/data/EdgeRTranscriptEstCounts.csv")
df$X <- NULL

#################### convert Ensembl transcript IDs to Ensembl gene IDS ######################
# convert to Ensembl Transcript ID to Ensembl Gene ID using biomaRt
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) # if unresponsive run: httr::set_config(httr::config(ssl_verifypeer = FALSE))
genes <- getBM(
  attributes=c("ensembl_gene_id", "ensembl_transcript_id"),
  mart=mart,
  useCache = FALSE)

# find the max count for different isoforms of the same gene
colnames(genes)[colnames(genes) == "ensembl_transcript_id"] <- "target_id"
df <- merge(df, genes, by="target_id")
df$target_id <- NULL
df <- df %>% group_by(df$ensembl_gene_id) %>% summarise(across(where(is.numeric), max))
df <- as.data.frame(df)
rownames(df) <- df$"df$ensembl_gene_id"
df$"df$ensembl_gene_id" <- NULL

# write.csv(df, "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN/results/data/EdgeRGeneEstCountsMax.csv")
df <- read.csv("/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN/results/data/EdgeRGeneEstCountsMax.csv")
rownames(df) <- df$X
df$X <- NULL

# create DGEList object
group <- factor(metaData$condition)
y <- DGEList(counts=df,group=group)

# create design matrix
batch <- metaData$LANE # test batch as covariate
# design <- model.matrix(~group)
design <- model.matrix(~group+batch)

# filtering some of the samples using default settings; read the documentation
# https://rdrr.io/bioc/edgeR/src/R/filterByExpr.R for source code
# below article says that from their results 14177 genes are typically expressed in kidneys
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2781110/
keep <- filterByExpr(y, design=design, group=group)
y <- y[keep,,keep.lib.sizes=FALSE]

# calculates effective library size using TMM (Trimmed Mean of M values)
y <- calcNormFactors(y)

# generate negative binomial likelihhod model using weighted empircal bayes
# calculates the likelihood by conditioning on the adjusted counts for each gene
y <- estimateDisp(y,design)

# LRT
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
DE_All <- as.data.frame(topTags(lrt, n=nrow(lrt$table)))
write.csv(DE_All, "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN/results/data/EdgeR_RNA_all_genes.csv")

# max isofrom transcript to gene counts
DE_All <- DE_All[DE_All$PValue < 3.53e-6,] # 3.53e-6 is equal to 0.05/14149
DE_sig <- DE_All[(DE_All$logFC > 2) | (DE_All$logFC < -2),] # filters to 478 genes with 1.5 and 185 with 2

write.csv(DE_sig, "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN/results/data/EdgeR_RNA_sig_genes.csv")

######################## convert Kallisto outputs using package ############################
# converts kallisto outputs to a table of counts (rows are transcript ids and columns are samples)
# makes same output as code at top of this file

# library(tximport)
# library(rhdf5)
# baseDir <- "/home/tnagano/projects/def-ccastel/SharedResources/TFAM_KO/FASTQ"

# metaData <- read.csv(paste(baseDir, "targetsSleuthRD.csv", sep="/")) # use to label samples
# folders <- paste0("kallisto_output", 1:12) # specify folders with kallisto outputs
# files <- file.path(baseDir, folders, "abundance.h5") # create a vector with file paths
# names(files) <- metaData$sample # label files with sample names
# txi <- tximport(files, type = "kallisto", txOut = TRUE) 

# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) # if unresponsive run: httr::set_config(httr::config(ssl_verifypeer = FALSE))
# genes <- getBM(
#   attributes=c("ensembl_transcript_id", "ensembl_gene_id"),
#   mart=mart,
#   useCache = FALSE)

# txi1 <- summarizeToGene(txi, genes)

######################## Old sleuth code ##############################
# setwd("results/data")
# # df <- read.csv("SleuthWaldTestResults.csv", header=T) # wald test results
# df <- read.csv("SleuthAllGenesAnnotatedRNASeqResultsGeneWise_cleaned.csv", header=T) # Likelihood test results
# names(df)[names(df) == "target_id"] <- "ensembl_gene_id"

# # use biomaRt to find mapping info between ensembl and entrez gene IDS
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) # if unresponsive run: httr::set_config(httr::config(ssl_verifypeer = FALSE))
# genes <- getBM(
#   attributes=c("ensembl_gene_id", "chromosome_name", "start_position"),
#   mart=mart,
#   useCache = FALSE)

# # merge data with entrez gene ids
# df_anno <- merge(df, genes, by = "ensembl_gene_id")
# df_sig <- filter(df_anno, pval < 3.59e-6) # 3.590149e-06
# df_sig <- df_sig %>% arrange(pval)