# Old file no longer used
# if need to install packages

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("minfi")

library(data.table)
library(dplyr)
library(minfi)
library(RColorBrewer)
library(limma)

# Update baseDir to proper location
baseDir <- "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN"

# set working directory to be folder containing root
setwd(baseDir)
baseDirRun1 <- paste(baseDir, "/data/run1", sep="")
baseDirRun2 <- paste(baseDir, "/data/run2", sep="")

# read in csv files into data.table objects with meta information
targetsRun1 <- read.metharray.sheet(baseDirRun1)
targetsRun2 <- read.metharray.sheet(baseDirRun2)

# remove chem samples from targetsRun2
targetsRun2 <- filter(targetsRun2, Group != "Chem")

# load in the raw intensity data from .idat files using targets
rgSetRun1 <- read.metharray.exp(targets = targetsRun1)
rgSetRun2 <- read.metharray.exp(targets = targetsRun2)

# calculate the detection p-values for each CpG in the sample
# used later for removing poor performing probes
detPRun1 <- detectionP(rgSetRun1)
detPRun2 <- detectionP(rgSetRun2)

# output the quality control reports
# Decode equivalent to Sample_Name and Status equal to Group between runs
#run1
qcReport(rgSetRun1, sampNames=targetsRun1$Sample_Name, sampGroup = targetsRun1$Group, pdf="./results/qcreportRun1.pdf")
#run2
qcReport(rgSetRun2, sampNames=targetsRun2$Decode, sampGroup = targetsRun2$Group, pdf="./results/qcreportRun2.pdf")

# preprocess using preprocessQuantile
mSetSqRun1 <- preprocessQuantile(rgSetRun1, quantileNormalize = TRUE, stratified = TRUE) 
mSetSqRun2 <- preprocessQuantile(rgSetRun2, quantileNormalize = TRUE, stratified = TRUE)

# visualize Run1 normalization
pdf("./results/Run1_norm.pdf")
par(mfrow=c(1,2))
densityPlot(getBeta(mSetSqRun1), sampGroups=targetsRun1$Group, main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targetsRun1$Group)), text.col=brewer.pal(8,"Dark2"))
densityPlot(rgSetRun1, sampGroups=targetsRun1$Group, main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targetsRun1$Group)), text.col=brewer.pal(8,"Dark2"))
dev.off()

# visualize Run2 normalization
pdf("./results/Run2_norm.pdf")
par(mfrow=c(1,2))
densityPlot(getBeta(mSetSqRun2), sampGroups=targetsRun2$Group, main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targetsRun2$Group)), text.col=brewer.pal(8,"Dark2")
densityPlot(rgSetRun2, sampGroups=targetsRun2$Group, main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targetsRun2$Group)), text.col=brewer.pal(8,"Dark2"))
dev.off()

# Filtering out poor performing probes using detection p values for the samples
# returns a p-value for each probe by comparing methylated and unmethylated signal to negative control based on normal distribution
# documentation recommends to remove all probes with det p-value < 0.01

# reorder rows in detection p-value to match with the probe names from the methylation data
detPRun1 <- detPRun1[match(featureNames(mSetSqRun1),rownames(detPRun1)),] 
detPRun2 <- detPRun2[match(featureNames(mSetSqRun2),rownames(detPRun2)),]

# find all probes that have failed in one or more samples based on p-value < 0.01
keepRun1 <- rowSums(detPRun1 < 0.01) == ncol(mSetSqRun1) 
keepRun2 <- rowSums(detPRun2 < 0.01) == ncol(mSetSqRun2) 
mSetSqFltRun1 <- mSetSqRun1[keepRun1,]
mSetSqFltRun2 <- mSetSqRun2[keepRun2,]
# table(keepRun1) # 15078 probes removed
# table(keepRun2) # 2786 probes removed

# remove Pidsley Cross reactive probes Pidsley et al. https://doi.org/10.1186/s13059-016-1066-1
pidsley_probes <- read.csv(file=paste(baseDir, "data/PidsleyCrossReactiveProbesEPIC.csv", sep="/"), stringsAsFactors=FALSE)

# get probes in runs which are in list of cross-reactive probes
pidsRun1 <- !(featureNames(mSetSqFltRun1) %in% pidsley_probes$TargetID)
pidsRun2 <- !(featureNames(mSetSqFltRun2) %in% pidsley_probes$TargetID)
mSetSqFltRun1 <- mSetSqFltRun1[pidsRun1,]
mSetSqFltRun2 <- mSetSqFltRun2[pidsRun2,]
# table(pidsRun1) # 42547 probes removed
# table(pidsRun2) # 43043 probes removed
# Dr. Castellani said aroudn 43254 should be removed; some dropped due to below signifance?

# remove probes affected by SNPs keep default 0 for maf because there should be no genetic variation in cell culture
# should Single-base-pair extension (SBE) SNPs be removed? no don't do this
# mSetSqFltRun1: Before: 808234 After: 780974
# mSetSqFltRun2: Before: 820030 After: 792231
mSetSqFltRun1 <- dropLociWithSnps(mSetSqFltRun1)
mSetSqFltRun2 <- dropLociWithSnps(mSetSqFltRun2)

############### Now that the data is cleaned separately; how to compare and combine data?

# output MDS plot of the methylation data to results for Run 1
# shows PC 1 separates the Control and KO lines nicely
pdf("./results/Run1_MDS.pdf")
par(mfrow=c(1,1))
plotMDS(getM(mSetSqFltRun1), top=1000, gene.selection="common", cex=0.8)
legend("right", legend=levels(factor(targetsRun1$Group)),
       cex=0.65, bg="white")
dev.off()

# output MDS plot of the methylation data to results for Run 2
# shows PC 1 separates the Control and KO lines nicely
pdf("./results/Run2_MDS.pdf")
par(mfrow=c(1,1))
plotMDS(getM(mSetSqFltRun2), top=1000, gene.selection="common", cex=0.8)
legend("right", legend=levels(factor(targetsRun2$Group)),
       cex=0.65, bg="white")
dev.off()

# need to get out beta or m values
# should we be using m-values or beta values
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
# says m values should be used for statistical analysis because consistent std across
# distribution for m values, but not for beta values
betaRun1 <- as.data.frame(getBeta(mSetSqFltRun1))
betaRun2 <- as.data.frame(getBeta(mSetSqFltRun2))
mValRun1 <- as.data.frame(getM(mSetSqFltRun1))
mValRun2 <- as.data.frame(getM(mSetSqFltRun2))

# rename betaRun1 columns with Sample_names
for(i in 1:6){
        old_col <- paste(targetsRun1$Slide[i], targetsRun1$Array[i], sep="_")
        names(betaRun1)[names(betaRun1) == old_col] <- targetsRun1$Sample_Name[i]
        names(mValRun1)[names(mValRun1) == old_col] <- targetsRun1$Sample_Name[i]
}

# rename betaRun2 columns with decode
for(i in 1:6){
        old_col <- paste(targetsRun2$Slide[i], targetsRun2$Array[i], sep="_")
        names(betaRun2)[names(betaRun2) == old_col] <- targetsRun2$Decode[i]
        names(mValRun2)[names(mValRun2) == old_col] <- targetsRun2$Decode[i]
}

############ get differentially methylated regions using linear model in limma
# Get knockout and normal groups
TFAM_KO_Run1 <- factor(targetsRun1$Group)
TFAM_KO_Run2 <- factor(targetsRun2$Group)

# create design matrix differentiating samples by KO or Normal
designRun1 <- model.matrix(~0+TFAM_KO_Run1, data=targetsRun1)
designRun2 <- model.matrix(~0+TFAM_KO_Run2, data=targetsRun2)

# use design matrix to fit m values 
fitRun1 <- lmFit(mValRun1, designRun1)
fitRun2 <- lmFit(mValRun2, designRun2)

# create contrast matrix to compare the differential methylation
contMatrixRun1 <- makeContrasts(TFAM_KO_Run1Knockout - TFAM_KO_Run1Normal, TFAM_KO_Run1Normal - TFAM_KO_Run1Knockout, levels=designRun1)
contMatrixRun2 <- makeContrasts(TFAM_KO_Run2Knockout - TFAM_KO_Run2Normal, TFAM_KO_Run2Normal - TFAM_KO_Run2Knockout, levels=designRun2)

# fit the contrast matrix to get CpGs with significant differential methylation
fitRun1_2 <- contrasts.fit(fitRun1, contMatrixRun1)
fitRun2_2 <- contrasts.fit(fitRun2, contMatrixRun2)

# not needed below?
fitRun1_2 <- eBayes(fitRun1_2)
fitRun2_2 <- eBayes(fitRun2_2)

############ find differentially methylated probes using dmpFinder
# get m values
mValRun1 <- getM(mSetSqFltRun1)
mValRun2 <- getM(mSetSqFltRun2)

# look more into dmpFinder to figure out best options
dmpRun1 <- dmpFinder(mValRun1, pheno=targetsRun1$Group, type="continuous")
dmpRun2 <- dmpFinder(mValRun2, pheno=targetsRun2$Group, type="continuous")

