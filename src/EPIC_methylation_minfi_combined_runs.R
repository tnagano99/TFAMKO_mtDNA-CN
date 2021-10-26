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
# should Single-base-pair extension (SBE) SNPs be removed? no don't do this
# mSetSqFlt: Before: 807862 After: 781118
mSetSqFlt <- dropLociWithSnps(mSetSqFlt, snps = c("CpG"))

# need to get out beta or m values
# should we be using m-values or beta values
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
# says m values should be used for statistical analysis because consistent std across
# distribution for m values, but not for beta values
beta <- as.data.frame(getBeta(mSetSqFlt))
# mVal <- as.data.frame(getM(mSetSqFlt))
mVal <- getM(mSetSqFlt)

# rename beta columns with decode
for(i in 1:24){
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

# write output files to results
write.csv(beta, "/project/6061316/tnagano/TFAMKO_mtDNA-CN/results/beta.csv")
write.csv(mVal, "/project/6061316/tnagano/TFAMKO_mtDNA-CN/results/mVal.csv")

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

# look more into dmpFinder to figure out best options
dmp <- dmpFinder(mVal, pheno=targets$Group, type="categorical")
write.csv(dmp, "/project/6061316/tnagano/TFAMKO_mtDNA-CN/results/dmp.csv")

########### find differentially methylated regions using DMRcate

# annotate the m-values to include info on genomic position etc.
myAnnotation <- cpg.annotate(object = mVal, datatype = "array", what = "M", analysis.type = "differential",
                             design = design, contrasts = TRUE, cont.matrix = contMatrix, fdr = 0.05,
                             coef = "TFAM_KO_Knockout - TFAM_KO_Normal", arraytype = "EPIC")

# using the annotation identify the diferentially methylated regions
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
rangesDMR <- extractRanges(DMRs)

# remove DMRs with regions with less than 10 CpGs
rangesDMR <- rangesDMR[(elementMetadata(rangesDMR)[, "no.cpgs"] >= 10)]

# export the DMR_ranges to csv file
df = as(rangesDMR, "data.frame")
write.table(df, "./results/DMRranges.csv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# to read in rangesDMR
# df2 = read.table("./results/DMRranges.csv", header = TRUE, sep = "\t", dec = ".")
# gr2 = makeGRangesFromDataFrame(df2, keep.extra.columns=TRUE)

############### GO/KEGG analysis on DMRs using MissMethyl

DMR_GO <- goregion(subset, all.cpg = rownames(mVal), collection = "GO", array.type = "EPIC")
DMR_GO <- DMR_GO %>% arrange(P.DE)

DMR_KEGG <- goregion(subset, all.cpg = rownames(mVal), collection = "KEGG", array.type = "EPIC")
DMR_KEGG <- DMR_KEGG %>% arrange(P.DE)
# note neuroactive ligand-receptor interaction is top hit for DMR using KEGG

########## Outputs to create visualizations #################################

# output the quality control reports
# Decode equivalent to Sample_Name and Status equal to Group between runs
qcReport(rgSet, sampNames=targets$Decode, sampGroup = targets$Group, pdf="./results/qcreport.pdf")

# visualize  normalization
pdf("./results/norm.pdf")
par(mfrow=c(1,2))
densityPlot(getBeta(mSetSq), sampGroups=targets$Group, main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Group)), text.col=brewer.pal(8,"Dark2"))
densityPlot(rgSet, sampGroups=targets$Group, main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Group)), text.col=brewer.pal(8,"Dark2"))
dev.off()

# output MDS plot of the methylation data to results for Run 1
# shows PC 1 separates the Control and KO lines nicely
pdf("./results/MDS.pdf")
par(mfrow=c(1,1))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", cex=0.8)
legend("right", legend=levels(factor(targets$Group)),
       cex=0.65, bg="white")
dev.off()

# plot the first DMR regions
# get the labels of the column based on the targets csv
groups <- c(Knockout="magenta", Normal="forestgreen")
cols <- groups[as.character(targets$Group)]

pdf("./results/DMR_1.pdf")
par(mfrow=c(1,1))
DMR.plot(ranges=rangesDMR, dmr=1, CpGs=mVal, what="M", arraytype="EPIC", phen.col=cols, genome = "hg19")
dev.off()