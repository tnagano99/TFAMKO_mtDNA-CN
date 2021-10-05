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
RGsetRun1 <- read.metharray.exp(targets = targetsRun1)
RGsetRun2 <- read.metharray.exp(targets = targetsRun2)