# if need to install packages

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("minfi")
BiocManager::install("missMethyl")

library(data.table)
library(dplyr)
library(minfi)
library(missMethyl)
# library(qqman)
# library(limma)

######################################## MissMethyl Probe Enrichment ########################
# Update baseDir to proper location
baseDir <- "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN"
setwd(baseDir)

# read in summary stats data
# update summary stats depending on which data to perform GO?KEGG
# summaryStats <- as.data.frame(read.csv("/home/tnagano/projects/def-ccastel/tnagano/EPIC/EPIClmerResults_DMP_ContinuousShrink.csv"))
summaryStats <- as.data.frame(read.csv(paste(baseDir, "/results/data/dmp_cont_shrink_10.csv", sep = "")))
summaryStats <- as.data.frame(read.csv(paste(baseDir, "/results/data/dmp_cont_10.csv", sep = "")))
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
write.csv(sigGO, paste(baseDir, "/results/data/GO_Cont_Shrink.csv", sep=""))

sigKEGG <- gometh(sig.cpg = sigCpGs, all.cpg = allCpGs, collection = "KEGG", array.type = "EPIC", plot.bias = TRUE, prior.prob = TRUE, anno = annotate)
sigKEGG <- sigKEGG %>% arrange(P.DE) # sort in order by P value
write.csv(sigKEGG, paste(baseDir, "/results/data/KEGG_Cont_Shrink.csv", sep=""))
