# if need to install packages

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("minfi")
BiocManager::install("DMRcate")

library(data.table)
library(dplyr)
library(minfi)
library(RColorBrewer)
library(limma)
library(DMRcate)
library(ggplot2)
library(missMethyl)

# Update baseDir to proper location
baseDir <- "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN"

# set working directory to be folder containing root
setwd(baseDir)
baseDirData <- paste(baseDir, "/data/run_combined", sep="")

# read in csv files into data.table objects with meta information
targets <- read.metharray.sheet(baseDirData)

# remove chem samples from targets
targets <- filter(targets, Group != "Chem")

# load in the raw intensity data from .idat files using targets
# force = TRUE because raw .idat files from run 2 have 826 less channels?
# checked for number of CpGs matched all same either combined or not
rgSet <- read.metharray.exp(targets = targets, force = TRUE)

# calculate the detection p-values for each CpG in the sample
# used later for removing poor performing probes
detP <- detectionP(rgSet)

# preprocess using preprocessFunnorm does within array normalization like preprocessSWAN
# 865859 CpG sites after preprocessing
mSetSq <- preprocessFunnorm(rgSet, bgCorr = TRUE, dyeCorr = TRUE, keepCN = TRUE, ratioConvert = TRUE)

# Filtering out poor performing probes using detection p values for the samples
# returns a p-value for each probe by comparing methylated and unmethylated signal to negative control based on normal distribution
# documentation recommends to remove all probes with det p-value < 0.01
# reorder rows in detection p-value to match with the probe names from the methylation data
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

# find all probes that have failed in one or more samples based on p-value < 0.01
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
mSetSqFlt <- mSetSq[keep,]

# table(keep) # 15465 probes removed

# remove Pidsley Cross reactive probes Pidsley et al. https://doi.org/10.1186/s13059-016-1066-1
pidsley_probes <- read.csv(file=paste(baseDir, "data/PidsleyCrossReactiveProbesEPIC.csv", sep="/"), stringsAsFactors=FALSE)

# get probes in runs which are in list of cross-reactive probes
pids <- !(featureNames(mSetSqFlt) %in% pidsley_probes$TargetID)
mSetSqFlt <- mSetSqFlt[pids,]
# table(pids) # 42532 probes removed

# remove probes affected by SNPs keep default 0 for maf because there should be no genetic variation in cell culture
# mSetSqFlt: Before: 807862 After: 781118 After w/o SBE: 780637 
mSetSqFlt <- dropLociWithSnps(mSetSqFlt, snps = c("CpG"))

# need to get out beta or m values
# should we be using m-values or beta values
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
# says m values should be used for statistical analysis because consistent std across
# distribution for m values, but not for beta values
beta <- as.data.frame(getBeta(mSetSqFlt))
mVal <- as.data.frame(getM(mSetSqFlt))

# rename beta columns with decode
for(i in 1:12){
        old_col <- paste(targets$Slide[i], targets$Array[i], sep="_")
        names(beta)[names(beta) == old_col] <- targets$Decode[i]
        names(mVal)[names(mVal) == old_col] <- targets$Decode[i]
		names(detP)[names(detP) == old_col] <- targets$Decode[i]
}

# mean filter check between runs and filter out rows based on mean differences on negative controls
beta["rowMeansRun1"] <- rowMeans(beta %>% select(contains("HEK293T") & ends_with("1")))
beta["rowMeansRun2"] <- rowMeans(beta %>% select(contains("HEK293T") & ends_with("2")))

# create scatter plot of rowMeansRun1 vs rowMeansRun2
threshold <- 0.15
beta['Threshold'] <- abs(beta$rowMeansRun1 - beta$rowMeansRun2) < threshold

# threshold value used and number of probes remaining before 781118 probes
# 0.2: 777150; 0.15: 769026; 0.1: 747004 
# c('cg26094004', 'cg26563141', 'cg08899667', 'cg21051031', 'cg14575356', 'cg23513930')
# cg26094004, cg21051031, cg23513930 removed in both
# cg26563141 removed when threshold is 0.1 not when 0.15
beta <- beta[abs(beta$rowMeansRun1 - beta$rowMeansRun2) < threshold, ]
beta <- beta[, !(names(beta) %in% c("rowMeansRun1", "rowMeansRun2", "Threshold"))]

# drop same rows in mVal for further analysis
mVal <- mVal[(row.names(mVal) %in% row.names(beta)), ]

# convert beta and m-values to matrices for further analysis
beta <- data.matrix(beta)
mVal <- data.matrix(mVal)

# write output files to results
# write.csv(beta, paste(baseDir, "/results/data/beta.csv", sep=""))
# write.csv(mVal, paste(baseDir, "/results/data/mVal.csv", sep=""))

# plot CpGs from run1 vs run2
# pdf("./results/plots/NC_Means_Scatter_0.15.pdf")
# par(mfrow=c(1,1))
# ggplot(beta, aes(x=rowMeansRun1, y=rowMeansRun2)) + geom_point(aes(color = factor(Threshold)))
# dev.off()

###################################### find differentially methylated probes using dmpFinder #############################################

# Find differentially methylated probes using categories of knockout or not
# not used below
# dmp <- dmpFinder(mVal, pheno=targets$Group, type="categorical", shrinkVar=FALSE)
# write.csv(dmp, paste(baseDir, "/results/data/dmp_shrink.csv", sep=""))

# Find differentially methylated probes using mtDNACN as continuous
dmp_cont <- dmpFinder(mVal, pheno=targets$mtDNACN, type="continuous", shrinkVar=TRUE)
write.csv(dmp_cont, paste(baseDir, "/results/data/dmp_cont.csv", sep=""))

###################################### find differentially methylated regions using DMRcate  #############################################

# Get knockout and normal groups and batches as factors
Type <- factor(targets$Group) # knockout or normal
Batch <- factor(targets$Slide) # run 1 or run 2 based on slide number
ID <- factor(substr(targets$Decode, 1, nchar(targets$Decode)-2)) # specify the 6 different samples (3 KO 3 NC)
design <- model.matrix(~Type+Batch) # specify design matrix Type and ID and linearly dependent so using Type and Batch as predictors explanation is below

# FIXED
# annotate the beta values to include info on genomic position etc.
# Coefficients not estimable: IDTFAM_KO_Clone_6
# https://support.bioconductor.org/p/39385/
# this is likely due to confounding between Type and ID because the ID gives information on KO or NC for each of the samples individually
# but Type gives the group of samples which are KO or NC
#### NOTE ### adjusted cpg.annotate function source code to use p-value instead of fdr (parameter still says fdr btw) 
# code is in code_snippets.R copy and run the code before running below line
myAnnotation <- cpg.annotate("array", beta, arraytype = "EPIC", analysis.type = "differential", design = design, coef = 2, what = "Beta", fdr = 1e-7)
# returns 9110 significant probes using pval < 1e-7

##### NOTE #####
# regions are determined by applying gaussian smoothing to the CpG-site test statistics using a given bandwidth set by lambda
# model the test statistic using chi-square distribution
# apply p-value adjustment using FDR to find significant CpG-sites
# group nearby CpG sites using lambda
# In this case our regions are determined by significant CpG sites within lambda (1000 in this case) base pairs of each other
# and there needs to be a minimum of 10 CpG sites to be considered a region based on our criteria 
# Stouffer: Stouffer summary transform of the individual CpG FDRs.
# HMFDR: Harmonic mean of the individual CpG FDRs.
# Fisher: Fisher combined probability transform of the individual CpG FDRs.
# using the annotation identify the diferentially methylated regions
# automatically uses pcutoff of 0.01 based on fdr set in cpg.annotate()
# use betacutoff of 0.05; removes CpG where beta shift is less than 0.05

# for (i in seq(5, 30, 5)){
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2, min.cpgs = 10, betacutoff = 0.05) # # fdr is 1 using same cutoffs get 4259 regions
# DMRs1 <- dmrcate(myAnnotation1, lambda=1000, C=2) # 106722 DMRs if fdr is 1 without cutoffs
rangesDMR <- extractRanges(DMRs, genome = "hg19") 
# removes ranges with mean methylation between groups is less than or equal to 0.01
rangesDMR <- rangesDMR[rangesDMR$meandiff > 0.05]
# rangesDMR <- rangesDMR[rangesDMR$Fisher < 1.17e-5]
rangesDMR <- as.data.frame(rangesDMR)

# filter out regions falling below bonferonni correction 0.05 / 4259 = 1.17e-5
rangesDMR <- filter(rangesDMR, Fisher < 1.17e-5) # 1.17e-5

write.csv(rangesDMR, paste(baseDir, "/results/data/DMRS_anno.csv", sep=""))
# write.csv(rangesDMR, paste(baseDir, "/results/data/DMRcate_", i, ".csv", sep=""))
# rangesDMR <- as.data.frame(rangesDMR)
# mean(rangesDMR[,"no.cpgs"])

# use goregion to find differentially methylated regions
DMR_sigGO <- goregion(rangesDMR, all.cpg = rownames(mVal), collection = "GO", array.type = "EPIC", plot.bias=TRUE)
DMR_sigGO <- DMR_sigGO %>% arrange(P.DE)
write.csv(DMR_sigGO, paste(baseDir, "/results/data/GO_DMRcate_pval.csv", sep=""))

DMR_sigKEGG <- goregion(rangesDMR, all.cpg = rownames(mVal), collection = "KEGG", array.type = "EPIC", plot.bias=TRUE)
DMR_sigKEGG <- DMR_sigKEGG %>% arrange(P.DE)
# note neuroactive ligand-receptor interaction is top hit for DMR using KEGG
write.csv(DMR_sigKEGG, paste(baseDir, "/results/data/KEGG_DMRcate_pval.csv", sep=""))
# }

# Convert to data frame for Manhattan plot
rangesDMR <- as.data.frame(rangesDMR)

# create Manhattan plot
rangesDMR$chr <- gsub("chr","",rangesDMR$seqnames)
rangesDMR<-subset(rangesDMR, (rangesDMR$chr!="0"))
rangesDMR<-subset(rangesDMR, (rangesDMR$chr!="X"))
rangesDMR<-subset(rangesDMR, (rangesDMR$chr!="Y"))
rangesDMR$chr <- as.numeric(rangesDMR$chr)

colnames(rangesDMR)[colnames(rangesDMR) == "HMFDR"] <- "P.Value"
manhattan(DMS=rangesDMR, filename="DMP", sig=4.61) # there are 2082 ranges, 0.05/2082 ~= 10^-4.62

# export the DMR_ranges to csv file
df = as(rangesDMR, "data.frame")
write.csv(df, paste(baseDir, "/results/data/DMRranges.csv", sep=""))
# write.table(df, "./results/DMRranges.csv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# to read in rangesDMR
# df2 = read.table("./results/data/DMRranges.csv", header = TRUE, sep = "\t", dec = ".")
# gr2 = makeGRangesFromDataFrame(df2, keep.extra.columns=TRUE)


############################################# Outputs to create visualizations ###################################################

# output the quality control reports
# Decode equivalent to Sample_Name and Status equal to Group between runs
qcReport(rgSet, sampNames=targets$Decode, sampGroup = targets$Group, pdf="./results/plots/qcreport.pdf")

# visualize  normalization
pdf("./results/plots/norm.pdf")
par(mfrow=c(1,2))
densityPlot(getBeta(mSetSq), sampGroups=targets$Group, main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Group)), text.col=brewer.pal(8,"Dark2"))
densityPlot(rgSet, sampGroups=targets$Group, main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Group)), text.col=brewer.pal(8,"Dark2"))
dev.off()

# output MDS plot of the methylation data to results for Run 1
# shows PC 1 separates the Control and KO lines nicely
pdf("./results/plots/MDS.pdf")
par(mfrow=c(1,1))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", cex=0.8)
legend("right", legend=levels(factor(targets$Group)),
       cex=0.65, bg="white")
dev.off()

# plot the first DMR regions
# get the labels of the column based on the targets csv
groups <- c(Knockout="magenta", Normal="forestgreen")
cols <- groups[as.character(targets$Group)]

for (i in 1:10) {
	pdf(paste0("./results/plots/DMR_", i, "_10e-4.pdf"))
	par(mfrow=c(1,1))
	DMR.plot(ranges=rangesDMR, dmr=i, CpGs=beta, what="Beta", arraytype="EPIC", phen.col=cols, genome = "hg19")
	dev.off()
}

# plot CpG beta values against mtDNA-CN
# using top CpGs from dmpFinder analysis using mtDNA-CN as continuous variable
CpGs <- rownames(dmp_cont)[1:20] 
CpG <- CpGs[20]
CpG_beta <- as.vector(beta[CpG,])
mtDNACN <- targets$mtDNACN
df <- data.frame(CpG_beta, mtDNACN)
Type <- factor(targets$Group)
# ID <- substr(targets$Decode, 1, nchar(targets$Decode)-2)
values <- coef(lm(formula = CpG_beta ~ mtDNACN, data = df))

pdf(paste0("./results/plots/mtDNACN-", CpG, ".pdf"))
par(mfrow=c(1,1))
ggplot(df, aes(x=mtDNACN, y=CpG_beta)) + geom_point(aes(color = Type)) + geom_abline(intercept = values[1], slope = values[2])
dev.off()

# create qqplot using qqman for probe analysis
library(qqman)

# summaryStats <- as.data.frame(read.csv(paste(baseDir, "/results/data/dmp_cont_shrink.csv", sep = "")))
# summaryStats <- as.data.frame(read.csv(paste(baseDir, "/results/data/dmp_cont.csv", sep = "")))
# summaryStats <- as.data.frame(read.csv(paste(baseDir, "/results/data/dmp_shrink.csv", sep = "")))
summaryStats<- as.data.frame(read.csv(paste(baseDir, "/results/data/Linear_Mixed_Model_lmerResults_Categorical.csv", sep = "")))
# summaryStats <- as.data.frame(read.csv(paste(baseDir, "/results/data/SleuthAllGenesAnnotatedRNASeqResultsGeneWise_cleaned.csv", sep = "")))
# summaryStats <- as.data.frame(read.csv("Linear_Mixed_Model_lmerResults.csv"))

pdf(paste0("./results/plots/QQ_LMM_Cat.pdf"))
par(mfrow=c(1,1))
qq(summaryStats$pval)
dev.off()

# create qqplot using qqman for region analysis

pdf(paste0("./results/plots/QQ_DMR.pdf"))
par(mfrow=c(1,1))
qq(rangesDMR$HMFDR)
dev.off()

