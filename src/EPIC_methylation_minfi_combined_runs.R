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

# output the quality control reports
# Decode equivalent to Sample_Name and Status equal to Group between runs
qcReport(rgSet, sampNames=targets$Decode, sampGroup = targets$Group, pdf="./results/qcreport.pdf")

# preprocess using preprocessFunnorm does within array normalization like preprocessSWAN
# 865859 CpG sites after preprocessing
mSetSq <- preprocessFunnorm(rgSet, bgCorr = TRUE, dyeCorr = TRUE, keepCN = TRUE, ratioConvert = TRUE)

# visualize  normalization
pdf("./results/norm.pdf")
par(mfrow=c(1,2))
densityPlot(getBeta(mSetSq), sampGroups=targets$Group, main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Group)), text.col=brewer.pal(8,"Dark2")
densityPlot(rgSet, sampGroups=targets$Group, main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Group)), text.col=brewer.pal(8,"Dark2"))
dev.off()

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
# Dr. Castellani said aroudn 43254 should be removed; some dropped due to below signifance?

# remove probes affected by SNPs keep default 0 for maf because there should be no genetic variation in cell culture
# should Single-base-pair extension (SBE) SNPs be removed? no don't do this
# mSetSqFlt: Before: 807862 After: 781118
mSetSqFlt <- dropLociWithSnps(mSetSqFlt, snps = c("CpG"))

############### Now that the data is cleaned separately; how to compare and combine data?

# output MDS plot of the methylation data to results for Run 1
# shows PC 1 separates the Control and KO lines nicely
pdf("./results/MDS.pdf")
par(mfrow=c(1,1))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", cex=0.8)
legend("right", legend=levels(factor(targets$Group)),
       cex=0.65, bg="white")
dev.off()

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
}

# mean filter check between runs and filter out rows based on mean differences on negative controls
beta["rowMeansRun1"] <- rowMeans(beta %>% select(contains("HEK293T") & ends_with("1")))
beta["rowMeansRun2"] <- rowMeans(beta %>% select(contains("HEK293T") & ends_with("2")))

# threshold value used and number of probes remaining before 781118 probes
# 0.2: 777150 0.1: 747004 
threshold <- 0.1
beta <- beta[abs(beta$rowMeansRun1 - beta$rowMeansRun2) < threshold, ]
beta <- beta[, !(names(beta) %in% c("rowMeansRun1", "rowMeansRun2"))]

# drop same rows in mVal for further analysis
mVal <- mVal[(row.names(mVal) %in% row.names(beta)), ]

############ get differentially methylated regions using linear model in limma using m-values
# Get knockout and normal groups
TFAM_KO_ <- factor(targets$Group)

# create design matrix differentiating samples by KO or Normal
design <- model.matrix(~0+TFAM_KO_, data=targets)

# use design matrix to fit m values 
fit <- lmFit(mVal, design)

# create contrast matrix to compare the differential methylation
contMatrix <- makeContrasts(TFAM_KO_Knockout - TFAM_KO_Normal, TFAM_KO_Normal - TFAM_KO_Knockout, levels=design)

# fit the contrast matrix to get CpGs with significant differential methylation
fit_2 <- contrasts.fit(fit, contMatrix)

# view the amount of signifcant CpGs
# summary(decideTests(fit_2))
############ find differentially methylated probes using dmpFinder
# get m values
mVal <- getM(mSetSqFlt)

# look more into dmpFinder to figure out best options
dmp <- dmpFinder(mVal, pheno=targets$Group, type="continuous")

