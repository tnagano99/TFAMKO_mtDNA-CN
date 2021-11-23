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
# mSetSqFlt: Before: 807862 After: 781118
mSetSqFlt <- dropLociWithSnps(mSetSqFlt, snps = c("CpG"))

# need to get out beta or m values
# should we be using m-values or beta values
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
# says m values should be used for statistical analysis because consistent std across
# distribution for m values, but not for beta values
beta <- as.data.frame(getBeta(mSetSqFlt))
mVal <- as.data.frame(getM(mSetSqFlt))
detP <- as.data.frame(detP)

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

pdf("./results/plots/NC_Means_Scatter_0.15.pdf")
par(mfrow=c(1,1))
ggplot(beta, aes(x=rowMeansRun1, y=rowMeansRun2)) + geom_point(aes(color = factor(Threshold)))
dev.off()

# threshold value used and number of probes remaining before 781118 probes
# 0.2: 777150; 0.15: 769026; 0.1: 747004 
beta <- beta[abs(beta$rowMeansRun1 - beta$rowMeansRun2) < threshold, ]
beta <- beta[, !(names(beta) %in% c("rowMeansRun1", "rowMeansRun2", "Threshold"))]

# drop same rows in mVal for further analysis
mVal <- mVal[(row.names(mVal) %in% row.names(beta)), ]

# drop detP rows that are not in remaining beta for local FEM analysis
detP <- detP[(row.names(detP) %in% row.names(beta)), ]

# write output files to results
write.csv(beta, paste(baseDir, "/results/data/beta.csv", sep=""))
write.csv(mVal, paste(baseDir, "/results/data/mVal.csv", sep=""))
write.csv(detP, paste(baseDir, "/results/data/detP.csv", sep=""))

# convert beta and m-values to matrices for further analysis
beta <- data.matrix(beta)
mVal <- data.matrix(mVal)

###################################### find differentially methylated probes using a linear mixed model ##################################

# code to run a linear model for each probe
Batch <- targets$Slide # run 1 or run 2 based on slide number
ID <- substr(targets$Decode, 1, nchar(targets$Decode)-2) # specify the 6 different samples (3 KO 3 NC)
mtDNACN <- targets$mtDNACN
EPICTFAMKO <- beta

library(lme4)
## Define the lm function
runlme <- function(thisdat) {
   lme1 <- eval(parse(text=expression));    
   ##Get the summary of the model
   smodel = summary(lme1);
   return(smodel)
}

varnames <- c("mtDNACN")
Bmat <- SEmat <- Tmat <- matrix(NA,nrow=nrow(EPICTFAMKO),ncol=1)
rownames(Bmat) <- rownames(SEmat) <- rownames(Tmat) <- rownames(EPICTFAMKO)
colnames(Bmat) <- paste("Estimate",varnames,sep=".")
colnames(SEmat) <- paste("Std.Error",varnames,sep=".")
colnames(Tmat) <- paste("t-value",varnames,sep=".")

for (i in 1:nrow(EPICTFAMKO)) { 
	if (i %% 10000 == 0) {
		cat(paste("On probe ",i,"\n",sep=""))
	} #outputs every 10000 probes just to check in

	thisExpr <- as.numeric(EPICTFAMKO[i,])
	expression <- "lmer(thisExpr~mtDNACN + (1|Batch)+(1|ID), na.action=na.exclude, control = lmerControl(calc.derivs = FALSE), REML=FALSE)"

	designmatrix <- data.frame(thisExpr, mtDNACN, Batch, ID)
	lme1.out <- try(runlme(designmatrix),silent=F);

	if (substr(lme1.out[1],1,5)!="Error") {
		tabOut <- lme1.out$coefficients
		Bmat[i,] <- tabOut[2,"Estimate"]
		SEmat[i,] <- tabOut[2,"Std. Error"]
		Tmat[i,] <- tabOut[2,"t value"]
	} else {
		cat('Error in LME of Probe',rownames(EPICTFAMKO)[i],"id",'\n')
		cat('Setting P-value=NA,Beta value=NA, and =NA\n');
		Bmat[i,] <- SEmat[i,] <- Tmat[i,] <- NA;
	}
}

warnings()

FinalResults <- cbind(Bmat, SEmat, Tmat)
FinalResults <- as.data.frame(FinalResults)
zscores <- FinalResults[,3]
pvalue <- pchisq(zscores**2,1,lower.tail=F)
min(pvalue)
FinalResults2 <- cbind(FinalResults,pvalue)

write.csv(FinalResults2, paste(baseDir, "results/Linear_Mixed_Model_lmerResults.csv", sep = "/"), quote=F)

###################################### find differentially methylated probes using dmpFinder #############################################

# Find differentially methylated probes using categories of knockout or not
dmp <- dmpFinder(mVal, pheno=targets$Group, type="categorical", shrinkVar=TRUE)
write.csv(dmp, paste(baseDir, "/results/data/dmp.csv", sep=""))

# Find differentially methylated probes using mtDNACN as continuous
dmp_cont <- dmpFinder(mVal, pheno=targets$mtDNACN, type="continuous", shrinkVar=TRUE)
write.csv(dmp_cont, paste(baseDir, "/results/data/dmp_cont.csv", sep=""))

###################################### find differentially methylated regions using DMRcate  #############################################

# Get knockout and normal groups and batches as factors
Type <- factor(targets$Group) # knockout or normal
Batch <- factor(targets$Slide) # run 1 or run 2 based on slide number
ID <- factor(substr(targets$Decode, 1, nchar(targets$Decode)-2)) # specify the 6 different samples (3 KO 3 NC)

# create design matrix differentiating samples, KO or Normal and batch effects
design <- model.matrix(~Type+Batch+ID) #1 some column vectors are linearly dependent
design <- model.matrix(~Batch+ID) #2
design <- model.matrix(~Type+ID) #3 some column vector are linearly dependent
design <- model.matrix(~Type+Batch) #4
# write.csv(design, paste(baseDir, "/results/data/design.csv", sep=""))

# annotate the beta values to include info on genomic position etc.
# Coefficients not estimable: IDTFAM_KO_Clone_6
# https://support.bioconductor.org/p/39385/
# this is likely due to confounding between Type and ID because the ID gives information on KO or NC for each of the samples individually
# but Type gives the group of samples which are KO or NC
myAnnotation <- cpg.annotate("array", beta, arraytype = "EPIC", analysis.type = "differential", design = design, coef = 2, what = "Beta", fdr = 0.01)

# using the annotation identify the diferentially methylated regions
# automatically uses pcutoff of 0.01 based on fdr set in cpg.annotate()
# use betacutoff of 0.05; removes 
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2, min.cpgs = 10, betacutoff = 0.05) 
rangesDMR <- extractRanges(DMRs, genome = "hg19") # 1) 3843 2) 15 3) 2644 4) 2256
testRanges <- rangesDMR[rangesDMR$meandiff > 0.01] # 1) 3659 2) 12 3) 2480 4) 2082

# create Manhattan plot
rangesDMR$chr <- gsub("chr","",rangesDMR$seqnames)
rangesDMR<-subset(rangesDMR, (rangesDMR$seqnames!="0"))
rangesDMR<-subset(rangesDMR, (rangesDMR$seqnames!="X"))
rangesDMR<-subset(rangesDMR, (rangesDMR$seqnames!="Y"))
rangesDMR$seqnames <- as.numeric(rangesDMR$seqnames)

colnames(rangesDMR)[colnames(rangesDMR) == "seqnames"] <- "chr"
colnames(rangesDMR)[colnames(rangesDMR) == "HMFDR"] <- "P.Value"
colnames(rangesDMR)[colnames(rangesDMR) == "P.Value"] <- "HMFDR"
colnames(rangesDMR)[colnames(rangesDMR) == "Fisher"] <- "P.Value"
manhattan(DMS=rangesDMR, filename="DMRcate", sig=6)

# export the DMR_ranges to csv file
df = as(rangesDMR, "data.frame")
write.csv(df, paste(baseDir, "/results/data/DMRranges.csv", sep=""))
# write.table(df, "./results/DMRranges.csv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# to read in rangesDMR
# df2 = read.table("./results/data/DMRranges.csv", header = TRUE, sep = "\t", dec = ".")
# gr2 = makeGRangesFromDataFrame(df2, keep.extra.columns=TRUE)

# use goregion to find differentially methylated regions
DMR_sigGO <- goregion(rangesDMR, all.cpg = rownames(mVal), collection = "GO", array.type = "EPIC", plot.bias=TRUE)
DMR_sigGO <- DMR_sigGO %>% arrange(P.DE)
write.csv(DMR_sigGO, paste(baseDir, "/results/data/GO_DMRcate.csv", sep=""))

DMR_sigKEGG <- goregion(rangesDMR, all.cpg = rownames(mVal), collection = "KEGG", array.type = "EPIC", plot.bias=TRUE)
DMR_sigKEGG <- DMR_sigKEGG %>% arrange(P.DE)
# note neuroactive ligand-receptor interaction is top hit for DMR using KEGG
write.csv(DMR_sigKEGG, paste(baseDir, "/results/data/KEGG_DMRcate.csv", sep=""))

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

pdf("./results/plots/DMR_1.pdf")
par(mfrow=c(1,1))
DMR.plot(ranges=rangesDMR, dmr=1, CpGs=beta, what="Beta", arraytype="EPIC", phen.col=cols, genome = "hg19")
dev.off()

