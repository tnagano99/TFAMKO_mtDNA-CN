# if need to install packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("minfi")
# BiocManager::install("DMRcate")

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
# mSetSqFlt: Before: 807862 After: 781118
mSetSqFlt <- dropLociWithSnps(mSetSqFlt, snps = c("CpG"))

# need to get out beta or m values
# should we be using m-values or beta values
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
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
beta["rowMeansRun1"] <- rowMeans(beta[,c("HEK293T_NC-1_1", "HEK293T_NC-2_1", "HEK293T_NC-3_1")])
beta["rowMeansRun2"] <- rowMeans(beta[,c("HEK293T_NC-1_2", "HEK293T_NC-2_2", "HEK293T_NC-3_2")])

# set threshold to remove CpGs with greater than difference between NC means
threshold <- 0.15
beta['Threshold'] <- abs(beta$rowMeansRun1 - beta$rowMeansRun2) < threshold

# create scatter plot of rowMeansRun1 vs rowMeansRun2
# pdf("./results/plots/minfi/NC_Means_Scatter_0.15.pdf")
# par(mfrow=c(1,1))
# ggplot(beta, aes(x=rowMeansRun1, y=rowMeansRun2)) + geom_point(aes(color = factor(Threshold)))
# dev.off()

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

###################################### find differentially methylated probes using dmpFinder #############################################

# Find differentially methylated probes using mtDNACN as continuous
dmp_cont <- dmpFinder(mVal, pheno=targets$mtDNACN, type="continuous", shrinkVar=TRUE)
write.csv(dmp_cont, paste(baseDir, "/results/data/dmp_cont.csv", sep=""))

###################################### find differentially methylated regions using DMRcate  #############################################
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
# Get knockout and normal groups and batches as factors
##### END NOTE #####

Type <- factor(targets$Group) # knockout or normal
Batch <- factor(targets$Slide) # run 1 or run 2 based on slide number
# ID <- factor(substr(targets$Decode, 1, nchar(targets$Decode)-2)) # specify the 6 different samples (3 KO 3 NC)
design <- model.matrix(~Type+Batch) # specify design matrix Type and ID and linearly dependent so using Type and Batch as predictors explanation is below

#################################
############## NOTE #############
#################################
# adjusted cpg.annotate function source code to use p-value instead of fdr (parameter still says fdr)
# code is in code_snippets.R under header "DMRcate cpg.annotate" 
# copy and run the code before running below line
myAnnotation <- cpg.annotate("array", beta, arraytype = "EPIC", analysis.type = "differential", design = design, coef = 2, what = "Beta", fdr = 1e-7)
# returns 9110 significant probes using pval < 1e-7

DMRs <- dmrcate(myAnnotation, lambda=1000, C=2, min.cpgs = 10, betacutoff = 0.05) # # fdr is 1 using same cutoffs get 4259 regions; 106722 DMRs if fdr is 1 without cutoffs
rangesDMR <- extractRanges(DMRs, genome = "hg19")

# removes ranges with mean methylation between groups is less than or equal to 0.01
rangesDMR <- rangesDMR[rangesDMR$meandiff > 0.05]
rangesDMR <- rangesDMR[rangesDMR$Fisher < 1.17e-5] # filter out regions falling below bonferonni correction 0.05 / 4259 = 1.17e-5

rangesDMR <- as.data.frame(rangesDMR)
# write.csv(rangesDMR, paste(baseDir, "/results/data/DMRS_anno.csv", sep=""))

############################################# Outputs to create visualizations ###################################################

############# Probe analysis viz from minfi ################
# output the quality control reports
# Decode equivalent to Sample_Name and Status equal to Group between runs
qcReport(rgSet, sampNames=targets$Decode, sampGroup = targets$Group, pdf="./results/plots/minfi/qcreport.pdf")

# visualize normalization
pdf("./results/plots/minfi/norm.pdf")
par(mfrow=c(1,2))
densityPlot(getBeta(mSetSq), sampGroups=targets$Group, main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Group)), text.col=brewer.pal(8,"Dark2"))
densityPlot(rgSet, sampGroups=targets$Group, main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Group)), text.col=brewer.pal(8,"Dark2"))
dev.off()

# output MDS plot of the methylation data to results for Run 1
# shows PC 1 separates the Control and KO lines nicely
pdf("./results/plots/minfi/MDS.pdf")
par(mfrow=c(1,1))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", cex=0.8)
legend("right", legend=levels(factor(targets$Group)),
       cex=0.65, bg="white")
dev.off()

############## DMR analysis viz from DMRcate ########
# plot the first DMR regions
# get the labels of the column based on the targets csv
groups <- c(Knockout="magenta", Normal="forestgreen")
cols <- groups[as.character(targets$Group)]

for (i in 1:10) {
	pdf(paste0("./results/plots/DMRcate/DMR_", i, "_10e-4.pdf"))
	par(mfrow=c(1,1))
	DMR.plot(ranges=rangesDMR, dmr=i, CpGs=beta, what="Beta", arraytype="EPIC", phen.col=cols, genome = "hg19")
	dev.off()
}
