#ELMER, FUNCTIONAL ENRICHMENT and HOMER
#if (!requireNamespace("BiocManager", quietly=TRUE))
  #  install.packages("BiocManager")
#BiocManager::install("ELMER")
#BiocManager::install("GenomicInteractions")

library(ELMER)
library(MultiAssayExperiment)

# setwd("/dcs01/arking/arkinglab/active/projects/aric/epigenetics/ccastell/EPIC/Data/EPICRerun/")
setwd("../results/data")
#Methylation for trusted probes (mean <0.2 between plates)
#Setup for ttest method using 0.2 dataset (least stringent of the three options seems to show the most biology)

#Make a note of all the filtering you did here*****************
EPIC <- read.csv("beta.csv", header=T)
row.names(EPIC) <- EPIC$X
EPIC$X <- NULL

# EPIC <- EPIC[,which(colnames(EPIC) %in% c("TFAM_KO_Clone_11_2", "HEK293T_NC.1_2", "TFAM_KO_Clone_5_2", "HEK293T_NC.2_2", "TFAM_KO_Clone_6_2", "HEK293T_NC.3_2", "HEK293T_NC.1_1", "HEK293T_NC.2_1", "HEK293T_NC.3_1", "TFAM_KO_Clone_5_1",  "TFAM_KO_Clone_6_1",  "TFAM_KO_Clone_11_1",))]

#file from RNASeq

RNA <- read.csv("SleuthGeneNormalized_TPM.csv", header=T)
row.names(RNA) <- RNA$X
RNA$X <- NULL
colnames(RNA) <- c("TFAM_KO_Clone_11_1", "TFAM_KO_Clone_11_2", "TFAM_KO_Clone_5_1", "TFAM_KO_Clone_5_2", "TFAM_KO_Clone_6_1", "TFAM_KO_Clone_6_2", "HEK293T_NC.1_1", "HEK293T_NC.1_2", "HEK293T_NC.2_1", "HEK293T_NC.2_2", "HEK293T_NC.3_1", "HEK293T_NC.3_2")
RNA <- RNA[, colnames(EPIC)]

#Remove genes with TPM < 0.5 in >49% of samples
RNA$TPM05 <- rowSums(RNA<0.5)
RNA <- subset(RNA, RNA$TPM05 < 7)
RNA$TPM05 <- NULL


primary <- colnames(RNA)
GroupLabel <- c("Experiment", "Control", "Experiment", "Control", "Experiment", "Control", "Control", "Control", "Control", "Experiment", "Experiment", "Experiment")
Sample <- cbind(primary, GroupLabel)
Sample <- as.data.frame(Sample)
row.names(Sample) <- Sample$primary

data <- createMAE(exp = RNA, 
                  met = EPIC,
                  met.platform = "EPIC",
                  genome = "hg19",
                  save = FALSE,
                  TCGA = FALSE, 
                  colData=Sample,
                  met.na.cut=0.05
                  )
#supervised mode


sig.diff <- get.diff.meth(data = data, 
                          group.col = "GroupLabel",
                          group1 =  "Control",
                          group2 = "Experiment",
                          mode = "supervised",
                          minSubgroupFrac = 1, # if supervised mode set to 1
                          sig.dif = 0,
                          diff.dir = "both",
                          cores = 1,  
                          pvalue = 0.05,
                          save = TRUE
                          )

#adjusted P value by BH
#ELMER will search for probes differently methylated in group Control (n:6) compared to Experiment (n:6)
#ooo Arguments ooo
#o Number of probes: 728012
#o Beta value difference cut-off: 0
#o FDR cut-off: 0.05
#o Mode: unsupervised
#o % of samples per group in each comparison: 1
#o Min number of samples per group in each comparison: 5
#o Nb of samples group1 in each comparison: 6
#o Nb of samples group2 in each comparison: 6
#71441 results

sig.diff <- sig.diff[order(sig.diff$pvalue),]

nearGenes <- GetNearGenes(data = data, 
                         probes = sig.diff$probe, 
                         numFlankingGenes = 40) # 20 upstream and 20 dowstream genes

pairs <- get.pair(data = data,
                      group.col = "GroupLabel",
                      group1 =  "Control",
                      group2 = "Experiment",
                      nearGenes = nearGenes,
                      mode = "supervised",
                      minSubgroupFrac = 1,
                      permu.size = 100000, # Please set to 100000 to get significant results
                      raw.pvalue = 0.05,   
                      Pe = 0.001, # Please set to 0.001 to get significant results
                      filter.probes = FALSE, # See preAssociationProbeFiltering function
                      filter.percentage = 0.05,
                      filter.portion = 0.3,
                      cores = 1,
                      label = "ALL", diff.dir="both")

save.image("ELMER_TFAMKO_INTER.RData")

pairs2 <- subset(pairs, pairs$Raw.p < 0.01)

enriched.motif <- get.enriched.motif(data = data,
                                     probes = pairs2$Probe, 
                                     label = "ALL",
                                     min.incidence = 10,
                                     lower.OR = 1.1,
                                     save = TRUE)

TF <- get.TFs(data = data, 
              group.col = "GroupLabel",
              group1 =  "Experiment",
              group2 = "Control",
              mode = "supervised",
              enriched.motif = enriched.motif,
              cores = 1, 
              label = "ALL", diff.dir="both")

save.image("ELMER_TFAMKO_FINAL.RData")


write.csv(pairs, "getPair.ALL.all.pairs.statistic_20.csv", quote=F)

#########


# scatter.plot(data = data,
#        byPair = list(probe = c("cg00095138"), gene = c("ENSG00000145423")), 
#        category = "GroupLabel", save = TRUE, lm_line = TRUE)

# #need to again check if you see this in the other direction ... all 0's in TFAM and expression in controls

# #One probe and all nearby genes
# scatter.plot(data = data,
#              byProbe = list(probe = c("cg14557185"), numFlankingGenes = 20), 
#              category = "GroupLabel", 
#              lm = TRUE, # Draw linear regression curve
#              save = TRUE) 

# #schematic plot nearby genes
# schematic.plot(pair = pairs, 
#                data = data,
#                group.col = "GroupLabel",
#                byProbe = pairs$Probe[1],
#                save = TRUE)

# #schematic plot nearby probes together
# schematic.plot(pair = pairs, 
#                data = data,   
#                group.col = "GroupLabel", 
#                byGene = pairs$GeneID[1],
#                save = TRUE)

# #TF expression vs. average DNA methylation
# scatter.plot(data = data,
#              byTF = list(TF = c("GATA1","GATA2"),
#                          probe = enriched.motif[[names(enriched.motif)[2]]]), 
#              category = "GroupLabel",
#              save = TRUE, 
#              lm_line = TRUE)

# enrich2 <- read.csv("getMotif.ALL.motif.enrichment.csv", header=T)

# motif.enrichment.plot(motif.enrichment = enrich2, 
#                       significant = list(OR = 1.1,lowerOR = 1.1), 
#                       label = "ALL", 
#                       summary = TRUE,
#                       save = FALSE)  

# heatmapPairs(data = data, 
#              group.col = "GroupLabel",
#               group1 =  "Experiment",
#               group2 = "Control",
#              pairs = pairs,
#              filename =  NULL)
             
# load("getTF.ALL.TFs.with.motif.pvalue.rda")
# motif <- colnames(TF.meth.cor)[1]
# TF.rank.plot(motif.pvalue = TF.meth.cor, 
#              motif = motif,
#              save = TRUE) 
             
             
             
#              pairs.first <- pair[match(unique(pair$Symbol), pair$Symbol),]
        
