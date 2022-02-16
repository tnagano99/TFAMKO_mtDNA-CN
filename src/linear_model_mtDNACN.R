library(lme4)
library(data.table)
library(dplyr)
library(minfi)
library(RColorBrewer)
library(limma)
library(DMRcate)

# script to run linear mixed model on beta values

# Update baseDir to proper location
baseDir <- "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN"

# set working directory to be folder containing root
setwd(baseDir)

baseDirData <- paste(baseDir, "/data/run_combined", sep="")

# read in csv files into data.table objects with meta information
targets <- read.metharray.sheet(baseDirData)

# remove chem samples from targets
targets <- filter(targets, Group != "Chem")

# code to run a linear model for each probe
# setting the information
Batch <- targets$Slide # run 1 or run 2 based on slide number
ID <- substr(targets$Decode, 1, nchar(targets$Decode)-2) # specify the 6 different samples (3 KO 3 NC)
# mtDNACN <- targets$mtDNACN # using mtDNA-CN
mtDNACN <- targets$Group # using categorical Knockout/Control

# read in the stored beta values csv
# set the first column to be row labels and then drop the column
EPICTFAMKO <- read.csv(paste(baseDir, "/results/data/beta.csv", sep=""))
rownames(EPICTFAMKO) <- EPICTFAMKO[,1]
EPICTFAMKO <- EPICTFAMKO[,2:ncol(EPICTFAMKO)]

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

write.csv(FinalResults2, paste(baseDir, "results/data/Linear_Mixed_Model_lmerResults_Categorical.csv", sep = "/"), quote=F)

###################################### find differentially methylated probes using a linear mixed model ##################################

# # code to run a linear model for each probe
# Batch <- targets$Slide # run 1 or run 2 based on slide number
# ID <- substr(targets$Decode, 1, nchar(targets$Decode)-2) # specify the 6 different samples (3 KO 3 NC)
# mtDNACN <- targets$mtDNACN
# EPICTFAMKO <- beta

# library(lme4)
# ## Define the lm function
# runlme <- function(thisdat) {
#    lme1 <- eval(parse(text=expression));    
#    ##Get the summary of the model
#    smodel = summary(lme1);
#    return(smodel)
# }

# varnames <- c("mtDNACN")
# Bmat <- SEmat <- Tmat <- matrix(NA,nrow=nrow(EPICTFAMKO),ncol=1)
# rownames(Bmat) <- rownames(SEmat) <- rownames(Tmat) <- rownames(EPICTFAMKO)
# colnames(Bmat) <- paste("Estimate",varnames,sep=".")
# colnames(SEmat) <- paste("Std.Error",varnames,sep=".")
# colnames(Tmat) <- paste("t-value",varnames,sep=".")

# for (i in 1:nrow(EPICTFAMKO)) { 
# 	if (i %% 10000 == 0) {
# 		cat(paste("On probe ",i,"\n",sep=""))
# 	} #outputs every 10000 probes just to check in

# 	thisExpr <- as.numeric(EPICTFAMKO[i,])
# 	expression <- "lmer(thisExpr~mtDNACN + (1|Batch)+(1|ID), na.action=na.exclude, control = lmerControl(calc.derivs = FALSE), REML=FALSE)"

# 	designmatrix <- data.frame(thisExpr, mtDNACN, Batch, ID)
# 	lme1.out <- try(runlme(designmatrix),silent=F);

# 	if (substr(lme1.out[1],1,5)!="Error") {
# 		tabOut <- lme1.out$coefficients
# 		Bmat[i,] <- tabOut[2,"Estimate"]
# 		SEmat[i,] <- tabOut[2,"Std. Error"]
# 		Tmat[i,] <- tabOut[2,"t value"]
# 	} else {
# 		cat('Error in LME of Probe',rownames(EPICTFAMKO)[i],"id",'\n')
# 		cat('Setting P-value=NA,Beta value=NA, and =NA\n');
# 		Bmat[i,] <- SEmat[i,] <- Tmat[i,] <- NA;
# 	}
# }

# warnings()

# FinalResults <- cbind(Bmat, SEmat, Tmat)
# FinalResults <- as.data.frame(FinalResults)
# zscores <- FinalResults[,3]
# pvalue <- pchisq(zscores**2,1,lower.tail=F)
# min(pvalue)
# FinalResults2 <- cbind(FinalResults,pvalue)

# write.csv(FinalResults2, paste(baseDir, "results/Linear_Mixed_Model_lmerResults.csv", sep = "/"), quote=F)