library(data.table)
library(dplyr)
library(minfi)
library(missMethyl)
library(biomaRt)
library(ggplot)
# library(qqman)
# library(limma)

# read in EPIC array annotation data
annotate <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
baseDir <- "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN"
######################################################## Manhattan Plots Probe and Region #####################################################

# create manhattan plot
# swap lmerResults and pval/pvalue depending on file
# lmerResults <- as.data.frame(read.csv(paste(baseDir, "/results/data/Linear_Mixed_Model_lmerResults.csv", sep = ""))) 
lmerResults <- as.data.frame(read.csv(paste(baseDir, "/results/data/dmp_cont.csv", sep = "")))
merged <- merge(lmerResults, annotate, by.x="X", by.y="Name")

colnames(merged)[colnames(merged) == "Estimate.mtDNACN"] <- "estimate"
# colnames(merged)[colnames(merged) == "pvalue"] <- "P.Value" # for linear mixed model
colnames(merged)[colnames(merged) == "pval"] <- "P.Value" # for dmpfinder using mtDNA-CN as continuous
colnames(merged)[colnames(merged) == "pos"] <- "start"

manhattanraw(DMS=merged, filename="DMP", sig=7)

# Manhattan plot for RNA data Sleuth
lmerResults <- as.data.frame(read.csv(paste(baseDir, "/results/data/SleuthAllGenesAnnotatedRNASeqResultsGeneWise_cleaned.csv", sep = "")))
colnames(lmerResults)[colnames(lmerResults) == "pval"] <- "P.Value" # for dmpfinder using mtDNA-CN as continuous
manhattanraw(DMS=lmerResults, filename="RNA", sig=5.445)

# Manhattan plot for RNA data EdgeR
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) # if unresponsive run: httr::set_config(httr::config(ssl_verifypeer = FALSE))
genes <- getBM(
  attributes=c("ensembl_gene_id", "chromosome_name", "start_position"),
  mart=mart,
  useCache = FALSE)

lmerResults <- as.data.frame(read.csv(paste(baseDir, "/results/data/EdgeR_RNA_all_genes.csv", sep = "")))
colnames(lmerResults)[colnames(lmerResults) == "X"] <- "ensembl_gene_id"
merged <- merge(lmerResults, genes, by="ensembl_gene_id", all = FALSE)
colnames(merged)[colnames(merged) == "PValue"] <- "P.Value"
colnames(merged)[colnames(merged) == "chromosome_name"] <- "chr"
colnames(merged)[colnames(merged) == "start_position"] <- "start"
chr_filt <- as.character(1:22)
merged <- filter(merged, merged$chr %in% chr_filt)
manhattanraw(DMS=merged, filename="RNA_EdgeR", sig=5.452)

############ DMRranges Manhattan plot ######################

ranges <- read.csv(paste(baseDir, "/results/data/DMRS_anno.csv", sep=""))
ranges$chr <- gsub("chr","",ranges$seqnames)
ranges<-subset(ranges, (ranges$seqnames!="0"))
ranges<-subset(ranges, (ranges$seqnames!="X"))
ranges<-subset(ranges, (ranges$seqnames!="Y"))
ranges$seqnames <- as.numeric(ranges$seqnames)

colnames(ranges)[colnames(ranges) == "seqnames"] <- "chr"
colnames(ranges)[colnames(ranges) == "HMFDR"] <- "P.Value"
manhattan(DMS=ranges, filename="DMRcate", sig=4.61) # there are 2082 ranges, 0.05/2082 ~= 10^-4.62

############### plot individual CpG beta values against mtDNA-CN #################
# using top CpGs from dmpFinder analysis using mtDNA-CN as continuous variable
dmp_cont <- as.data.frame(read.csv("./results/data/dmp_cont.csv"))
rownames(dmp_cont) <- dmp_cont$X
dmp_cont$X <- NULL
CpGs <- rownames(dmp_cont)[1:20] 
CpG <- CpGs[1]
CpG_beta <- as.vector(beta[CpG,])
mtDNACN <- targets$mtDNACN
df <- data.frame(CpG_beta, mtDNACN)
Type <- factor(targets$Group)
# ID <- substr(targets$Decode, 1, nchar(targets$Decode)-2)
values <- coef(lm(formula = CpG_beta ~ mtDNACN, data = df))

pdf(paste0("./results/plots/mtDNACN-", CpG, ".pdf"))
par(mfrow=c(1,1))
ggplot(df, aes(x=mtDNACN, y=CpG_beta)) + geom_point(aes(color = Type, size = 10)) + geom_abline(intercept = values[1], slope = values[2]) + ylab(CpG)
dev.off()

############# create qqplot using qqman for probe/region analysis ###################
library(qqman)

summaryStats <- as.data.frame(read.csv(paste(baseDir, "/results/data/dmp_cont.csv", sep = "")))
# summaryStats <- as.data.frame(read.csv(paste(baseDir, "/results/data/SleuthAllGenesAnnotatedRNASeqResultsGeneWise_cleaned.csv", sep = "")))

pdf(paste0("./results/plots/QQ_LMM_Cat.pdf"))
par(mfrow=c(1,1))
qq(summaryStats$pval)
dev.off()

# create qqplot using qqman for region analysis
rangesDMR <- 
pdf(paste0("./results/plots/QQ_DMR.pdf"))
par(mfrow=c(1,1))
qq(rangesDMR$HMFDR)
dev.off()