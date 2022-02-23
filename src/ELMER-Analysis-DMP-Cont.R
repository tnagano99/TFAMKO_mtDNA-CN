library(ELMER)
library(MultiAssayExperiment)
library(data.table)
library(dplyr)

# setwd("/dcs01/arking/arkinglab/active/projects/aric/epigenetics/ccastell/EPIC/Data/EPICRerun/")
setwd("../results/data")

load("ELMER_TFAMKO_FINAL.RData")

sig_cpgs <- as.data.frame(read.csv("dmp_cont.csv"))
sig_cpgs <- filter(sig_cpgs, pval < 1e-7)

nearGenes <- GetNearGenes(data = data, 
                         probes = sig_cpgs$X, 
                         numFlankingGenes = 20)

pairs <- get.pair(data = data,
                      group.col = "GroupLabel",
                      group1 =  "Control",
                      group2 = "Experiment",
                      nearGenes = nearGenes,
                      mode = "unsupervised",
                      minSubgroupFrac = 0.4,
                      permu.size = 100000, # Please set to 100000 to get significant results
                      raw.pvalue = 0.05,   
                      Pe = 0.001, # Please set to 0.001 to get significant results
                      filter.probes = FALSE, # See preAssociationProbeFiltering function
                      filter.percentage = 0.2,
                      filter.portion = 0.3,
                      cores = 1,
                      label = "ALL_DMP_CONT_20_UNSUPERVISED_NO_FILTER", diff.dir="both")

save.image("ELMER_TFAMKO_INTER_DMP_CONT_20_UNSUPERVISED_NO_FILTER.RData")

pairs <- read.csv("getPair.ALL_DMP_CONT_20_UNSUPERVISED_NO_FILTER.all.pairs.statistic.csv")
# pairs1 <- read.csv("getPair.ALL_DMP_CONT_20_UNSUPERVISED.pairs.significant.csv")
pairs2 <- subset(pairs, pairs$Raw.p < 0.01) # 8895 pairs p < 0.005

enriched.motif <- get.enriched.motif(data = data,
                                     probes = pairs2$Probe, 
                                     label = "ALL_DMP_CONT_20_UNSUPERVISED_NO_FILTER",
                                     min.incidence = 10,
                                     lower.OR = 1.1,
                                     save = TRUE)

TF <- get.TFs(data = data, 
              group.col = "GroupLabel",
              group1 =  "Experiment",
              group2 = "Control",
              mode = "unsupervised",
              enriched.motif = enriched.motif,
              cores = 1, 
              label = "ELMER_TFAMKO_FINAL_DMP_CONT_20_UNSUPERVISED_NO_FILTER", diff.dir="both")


save.image("ALL_DMP_CONT_20_UNSUPERVISED_NO_FILTER.RData")