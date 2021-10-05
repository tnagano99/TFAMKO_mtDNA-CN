if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("minfi")

library(data.table)
library(dplyr)
library(minfi)

# set working directory to be folder containing run1 and run2 data
baseDirRun1 <- paste(getwd(), "/run1", sep="")
baseDirRun2 <- paste(getwd(), "/run2", sep="")

# read in csv files into data.table objects with meta information
targetsRun1 <- read.metharray.sheet(baseDirRun1)
targetsRun2 <- read.metharray.sheet(baseDirRun2)

# remove chem samples from targetsRun2
targetsRun2 <- filter(targetsRun2, Group != "Chem")

# load in the raw intensity data from .idat files using targets
rgSetRun1 <- read.metharray.exp(targets = targetsRun1)
rgSetRun2 <- read.metharray.exp(targets = targetsRun2)

# calculate the detection p-values for each CpG in the sample
detPRun1 <- detectionP(rgSetRun1)
detPRun2 <- detectionP(rgSetRun2)

# remove samples from runs with detection p-value less than 0.05
# no samples above threshold
keep1 <- colMeans(detPRun1) < 0.05
keep2 <- colMeans(detPRun2) < 0.05

# output the quality control reports
# Decode equivalent to Sample_Name and Status equal to Group between runs
#run1
qcReport(rgSetRun1, sampNames=targetsRun1$Sample_Name, sampGroup = targetsRun1$Status, pdf="qcreportRun1.pdf")
#run2
qcReport(rgSetRun2, sampNames=targetsRun2$Decode, sampGroup = targetsRun2$Group, pdf="qcreportRun2.pdf")

# preprocess using preprocessQuantile
mSetSqRun1 <- preprocessQuantile(rgSetRun1) # returning error "Error in if (ncol(df) >0) { : argument is of length zero"
mSetSqRun2 <- preprocessQuantile(rgSetRun2)

# visualize Run2 normalization
pdf("Run2.pdf")
par(mfrow=c(1,2))
densityPlot(getBeta(mSetSqRun2), sampGroups=targetsRun2$Group, main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targetsRun2$Group)))
densityPlot(rgSetRun2, sampGroups=targetsRun2$Group, main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targetsRun2$Group)))
dev.off()