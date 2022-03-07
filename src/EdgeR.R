library(edgeR)
library(dplyr)
library(data.table)
library(biomaRt)

# getting estimated counts from Kallisto output
baseDir <- "/home/tnagano/projects/def-ccastel/SharedResources/TFAM_KO/FASTQ"

metaData <- read.csv(paste(baseDir, "targetsSleuthRD.csv", sep="/"))
df <- read.csv(paste(baseDir, "kallisto_output1/abundance.tsv", sep="/"), sep = "\t")
df <- df[,c("target_id","est_counts")]
colnames(df)[colnames(df) == "est_counts"] <- metaData$sample[1]

for (i in 2:12){
    subDir <- paste0("kallisto_output", i)
    tmp <- read.csv(paste(baseDir, subDir,"abundance.tsv", sep="/"), sep = "\t")
    tmp <- tmp[,c("target_id","est_counts")]
    colnames(tmp)[colnames(tmp) == "est_counts"] <- metaData$sample[i]
    df <- merge(df, tmp, by = "target_id")
}

# write.csv(df, "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN/results/data/SleuthTranscriptEstCounts.csv")
df <- read.csv("/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN/results/data/SleuthTranscriptEstCounts.csv")
df$X <- NULL

# convert to Gene ID
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) # if unresponsive run: httr::set_config(httr::config(ssl_verifypeer = FALSE))
genes <- getBM(
  attributes=c("ensembl_gene_id", "ensembl_transcript_id"),
  mart=mart,
  useCache = FALSE)

colnames(genes)[colnames(genes) == "ensembl_transcript_id"] <- "target_id"
df <- merge(df, genes, by="target_id")
df$target_id <- NULL
df <- df %>% group_by(df$ensembl_gene_id) %>% summarise(across(where(is.numeric), sum))
df <- as.data.frame(df)
rownames(df) <- df$"df$ensembl_gene_id"
df$"df$ensembl_gene_id" <- NULL

# write.csv(df, "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN/results/data/SleuthGeneEstCounts.csv")

group <- factor(metaData$condition)
y <- DGEList(counts=df,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
# fit <- glmFit(y,design)
# lrt <- glmLRT(fit,coef=2)
# topTags(lrt)

fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
DE_All <- as.data.frame(topTags(qlf, n=nrow(qlf$table))) #3.44e-6 is equal to 0.05/14527
DE_sig <- filter(DE_All, DE_All$PValue < 3.44e-6) # 367 significant genes should there be a logFC filter

setwd("results/data")
# df <- read.csv("SleuthWaldTestResults.csv", header=T) # wald test results
df <- read.csv("SleuthAllGenesAnnotatedRNASeqResultsGeneWise_cleaned.csv", header=T) # Likelihood test results
names(df)[names(df) == "target_id"] <- "ensembl_gene_id"

# use biomaRt to find mapping info between ensembl and entrez gene IDS
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) # if unresponsive run: httr::set_config(httr::config(ssl_verifypeer = FALSE))
genes <- getBM(
  attributes=c("ensembl_gene_id", "chromosome_name", "start_position"),
  mart=mart,
  useCache = FALSE)

# merge data with entrez gene ids
df_anno <- merge(df, genes, by = "ensembl_gene_id")
df_sig <- filter(df_anno, pval < 3.59e-6) # 3.590149e-06