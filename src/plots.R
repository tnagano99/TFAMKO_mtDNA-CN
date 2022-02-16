library(data.table)
library(dplyr)
library(minfi)
library(missMethyl)
# library(qqman)
# library(limma)

# read in EPIC array annotation data
annotate <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

######################################################## Manhattan Plots Probe and Region #####################################################

# create manhattan plot
# swap lmerResults and pval/pvalue depending on file
# lmerResults <- as.data.frame(read.csv(paste(baseDir, "/results/data/Linear_Mixed_Model_lmerResults.csv", sep = ""))) 
lmerResults <- as.data.frame(read.csv(paste(baseDir, "/results/data/dmp_cont_10.csv", sep = ""))) 
merged <- merge(lmerResults, annotate, by.x="X", by.y="Name")

colnames(merged)[colnames(merged) == "Estimate.mtDNACN"] <- "estimate"
# colnames(merged)[colnames(merged) == "pvalue"] <- "P.Value" # for linear mixed model
colnames(merged)[colnames(merged) == "pval"] <- "P.Value" # for dmpfinder using mtDNA-CN as continuous
colnames(merged)[colnames(merged) == "pos"] <- "start"

manhattanraw(DMS=merged, filename="DMP_10", sig=7.3)

# Manhattan plot for RNA data
lmerResults <- as.data.frame(read.csv(paste(baseDir, "/results/data/SleuthAllGenesAnnotatedRNASeqResultsGeneWise_cleaned.csv", sep = "")))
colnames(lmerResults)[colnames(lmerResults) == "pval"] <- "P.Value" # for dmpfinder using mtDNA-CN as continuous
manhattanraw(DMS=lmerResults, filename="RNA", sig=7.3)

# gsub("chr","",merged$chr)
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