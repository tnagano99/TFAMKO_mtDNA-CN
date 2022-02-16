# if need to install packages

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("minfi")
BiocManager::install("missMethyl")

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(minfi)
library(missMethyl)
library(pathview)
library(biomaRt)
# library(qqman)
# library(limma)

######################################## MissMethyl Probe Enrichment ########################
# Update baseDir to proper location
baseDir <- "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN"
setwd(baseDir)

# read in summary stats data
# update summary stats depending on which data to perform GO?KEGG
# summaryStats <- as.data.frame(read.csv("/home/tnagano/projects/def-ccastel/tnagano/EPIC/EPIClmerResults_DMP_ContinuousShrink.csv"))
# summaryStats <- as.data.frame(read.csv(paste(baseDir, "/results/data/dmp_cont_shrink.csv", sep = "")))
summaryStats <- as.data.frame(read.csv(paste(baseDir, "/results/data/dmp_cont.csv", sep = "")))
# summaryStats <- as.data.frame(read.csv(paste(baseDir, "/results/data/dmp.csv", sep = "")))
# summaryStats<- as.data.frame(read.csv(paste(baseDir, "/results/data/Linear_Mixed_Model_lmerResults.csv", sep = "")))

# read in EPIC array annotation data
annotate <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# find signifcant CpGs (used p-value < 1e-7)
# swap below pval for pvalue if using Linear_Mixed_Model
df_sigCpGs <- filter(summaryStats, pval < 1e-7)
sigCpGs <- df_sigCpGs$X
allCpGs <- rownames(annotate)

# use gometh with prior probabilities to account for bias of CpG sites per gene
sigGO <- gometh(sig.cpg = sigCpGs, all.cpg = allCpGs, collection = "GO", array.type = "EPIC", plot.bias = TRUE, prior.prob = TRUE, anno = annotate)
sigGO <- sigGO %>% arrange(P.DE) # sort in order by P value
write.csv(sigGO, paste(baseDir, "/results/data/GO_Cont.csv", sep=""))

sigKEGG <- gometh(sig.cpg = sigCpGs, all.cpg = allCpGs, collection = "KEGG", array.type = "EPIC", plot.bias = TRUE, prior.prob = TRUE, anno = annotate)
sigKEGG <- sigKEGG %>% arrange(P.DE) # sort in order by P value
write.csv(sigKEGG, paste(baseDir, "/results/data/KEGG_Cont.csv", sep=""))

# output differentially methylated probes with annotation

merged <- merge(df_sigCpGs, annotate, by.x="X", by.y="Name")
merged <- as.data.frame(merged)
merged <- merged %>% arrange(pval)
merged <- merged[,c("X", "intercept", "beta", "t", "pval", "qval", "chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group", "GencodeCompV12_NAME", "GencodeCompV12_Accession", "GencodeCompV12_Group")]
write.csv(merged, paste(baseDir, "/results/data/DMPs_anno.csv", sep=""))

############ KEGG pathview

# differentially methylated probes
subset <- as.data.frame(merged[,c("X","GencodeCompV12_Accession")])
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

pathways <- pathview(gene.data = expanded$ensembl_gene_id, pathway.id = "04727", species = "hsa", gene.idtype="ENSEMBLTRANS") # 04080 Neuro 04727 GABA synapse

setwd("results/data")
# differentially expressed genes
df <- read.csv("SleuthAllGenesAnnotatedRNASeqResultsGeneWise_cleaned.csv", header=T) # Likelihood test results
names(df)[names(df) == "target_id"] <- "ensembl_gene_id"

df <- filter(df, pval < 3.59e-6)
pathways <- pathview(gene.data = df$ensembl_gene_id, pathway.id = "04727", species = "hsa", gene.idtype="ENSEMBL", out.suffix = "RNA")

# differentially methylated region
df <- read.csv("DMRS_anno.csv", header=T) # Likelihood test results
df <- separate_rows(df, overlapping.genes, sep=", ", convert = TRUE)

pathways <- pathview(gene.data = df$overlapping.genes, pathway.id = "04080", species = "hsa", gene.idtype="SYMBOL", out.suffix = "DMR")