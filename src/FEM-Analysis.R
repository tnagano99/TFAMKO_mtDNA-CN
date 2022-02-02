library(FEM)

# set working directory
setwd("C:/Users/tnaga/Documents/R_Scripts/Thesis")

# read in beta and detP values from minfi and read in Sleuth output
beta <- read.csv("beta.csv")
detP <- read.csv("detP.csv")
sleuth_results <- read.csv("SleuthWaldTestResults.csv")

# set row names of beta and detP to be CpG sites
row.names(beta) <- beta$X
beta$X <- NULL
row.names(detP) <- detP$X
detP$X <- NULL
sleuth_results$X <- NULL