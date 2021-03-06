library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(minfi)
library(missMethyl)
library(biomaRt)
library(GenomicRanges)
library("org.Hs.eg.db")

############################ Functional Enrichment DMP MissMethyl ########################
# Update baseDir to proper location
baseDir <- "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN"
setwd(baseDir)

# read in summary stats data
# update summary stats depending on which data to perform GO/KEGG
summaryStats <- as.data.frame(read.csv(paste(baseDir, "/results/data/dmp_cont.csv", sep = "")))

# read in EPIC array annotation data
annotate <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# find signifcant CpGs (used p-value < 1e-7)
# swap below pval for pvalue if using Linear_Mixed_Model
df_sigCpGs <- filter(summaryStats, pval < 1e-7)
sigCpGs <- df_sigCpGs$X
allCpGs <- summaryStats$X # use for DMR below

# use gometh with prior probabilities to account for bias of CpG sites per gene
sigGO <- gometh(sig.cpg = sigCpGs, all.cpg = allCpGs, collection = "GO", array.type = "EPIC", plot.bias = TRUE, prior.prob = TRUE, anno = annotate)
sigGO <- sigGO %>% arrange(P.DE) # sort in order by P value
write.csv(sigGO, paste(baseDir, "/results/data/GO_Cont.csv", sep=""))

sigKEGG <- gometh(sig.cpg = sigCpGs, all.cpg = allCpGs, collection = "KEGG", array.type = "EPIC", plot.bias = TRUE, prior.prob = TRUE, anno = annotate)
sigKEGG <- sigKEGG %>% arrange(P.DE) # sort in order by P value
write.csv(sigKEGG, paste(baseDir, "/results/data/KEGG_Cont.csv", sep=""))

# output differentially methylated probes with annotation
merged <- merge(summaryStats, annotate, by.x="X", by.y="Name")
merged <- as.data.frame(merged)
merged <- merged %>% arrange(pval)
merged <- merged[,c("X", "intercept", "beta", "t", "pval", "qval", "chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group", "GencodeCompV12_NAME", "GencodeCompV12_Accession", "GencodeCompV12_Group")]
write.csv(merged, paste(baseDir, "/results/data/DMPs_anno.csv", sep=""))

####################### Functional Enrichment DMR MissMethyl ##############################

# read in DMRs
rangesDMR <- as.data.frame(read.csv(paste(baseDir, "/results/data/DMRS_anno.csv", sep = "")))
rangesDMR$X <- NULL

# convert dataframe to granges object
rangesDMR <- makeGRangesFromDataFrame(rangesDMR, keep.extra.columns = TRUE)

# use goregion to find differentially methylated regions
DMR_sigGO <- goregion(rangesDMR, all.cpg = allCpGs, collection = "GO", array.type = "EPIC", plot.bias=TRUE)
DMR_sigGO <- DMR_sigGO %>% arrange(P.DE)
write.csv(DMR_sigGO, paste(baseDir, "/results/data/GO_DMRcate.csv", sep=""))

DMR_sigKEGG <- goregion(rangesDMR, all.cpg = allCpGs, collection = "KEGG", array.type = "EPIC", plot.bias=TRUE)
DMR_sigKEGG <- DMR_sigKEGG %>% arrange(P.DE)
write.csv(DMR_sigKEGG, paste(baseDir, "/results/data/KEGG_DMRcate.csv", sep=""))

