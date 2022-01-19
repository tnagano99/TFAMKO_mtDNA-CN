# if need to install packages

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("minfi")
BiocManager::install("missMethyl")

library(data.table)
library(dplyr)
library(minfi)
library(missMethyl)
library(qqman)
# library(limma)

# Update baseDir to proper location
baseDir <- "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN"
setwd(baseDir)

# read in summary stats data
# update summary stats depending on which data to perform GO?KEGG
# summaryStats <- as.data.frame(read.csv("/home/tnagano/projects/def-ccastel/tnagano/EPIC/EPIClmerResults_DMP_ContinuousShrink.csv"))
summaryStats <- as.data.frame(read.csv(paste(baseDir, "/results/data/dmp_cont.csv", sep = "")))
summaryStats <- as.data.frame(read.csv(paste(baseDir, "/results/data/dmp.csv", sep = "")))
summaryStats<- as.data.frame(read.csv(paste(baseDir, "/results/data/Linear_Mixed_Model_lmerResults.csv", sep = "")))

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

# create manhattan plot
# swap lmerResults and pval/pvalue depending on file
lmerResults <- as.data.frame(read.csv(paste(baseDir, "/results/data/Linear_Mixed_Model_lmerResults.csv", sep = "")))
# lmerResults <- as.data.frame(read.csv(paste(baseDir, "/results/data/dmp_cont.csv", sep = "")))
merged <- merge(lmerResults, annotate, by.x="X", by.y="Name")

colnames(merged)[colnames(merged) == "Estimate.mtDNACN"] <- "estimate"
colnames(merged)[colnames(merged) == "pvalue"] <- "P.Value" # for linear mixed model
# colnames(merged)[colnames(merged) == "pval"] <- "P.Value" # for dmpfinder using mtDNA-CN as continuous
colnames(merged)[colnames(merged) == "pos"] <- "start"

manhattan(DMS=merged, filename="lmerResults_LMM", sig=7.3)

# qqman manhattan plot
# merged$chr <- gsub("chr","",merged$chr)
# merged<-subset(merged, (merged$chr!="0"))
# merged<-subset(merged, (merged$start!="NA"))
# merged<-subset(merged, (merged$chr!="X"))
# merged<-subset(merged, (merged$chr!="Y"))
# merged$chr <- as.numeric(merged$chr)

# colnames(merged)[colnames(merged) == "chr"] <- "CHR"
# colnames(merged)[colnames(merged) == "start"] <- "BP"
# # colnames(merged)[colnames(merged) == "P.Value"] <- "P"
# colnames(merged)[colnames(merged) == "pval"] <- "P"
# colnames(merged)[colnames(merged) == "Probe_rs"] <- "SNP"

# jpeg("./results/plots/manhattan_qqman_cont.jpeg",res=400,width = 50, height = 20,units="cm")
# manhattan(merged, genomewideline = -log10(1e-7), suggestiveline = FALSE)
# dev.off()

# DMRranges Manhattan plot
ranges$chr <- gsub("chr","",ranges$seqnames)
ranges<-subset(ranges, (ranges$seqnames!="0"))
ranges<-subset(ranges, (ranges$seqnames!="X"))
ranges<-subset(ranges, (ranges$seqnames!="Y"))
ranges$seqnames <- as.numeric(ranges$seqnames)

colnames(ranges)[colnames(ranges) == "seqnames"] <- "chr"
colnames(ranges)[colnames(ranges) == "HMFDR"] <- "P.Value"
manhattan(DMS=ranges, filename="DMRcate", sig=7.3)

