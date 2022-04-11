#ELMER, FUNCTIONAL ENRICHMENT and HOMER
library(ELMER)
library(MultiAssayExperiment)
library(data.table)
library(dplyr)
library(biomaRt)
library(sesameData)
library(matrixStats)

# update to your working directory
baseDir <- "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN"
setwd(baseDir)

######################### ELMER Creating Data Input ##########################
# Read in beta values
EPIC <- read.csv("./results/data/beta.csv", header=T)
row.names(EPIC) <- EPIC$X
EPIC$X <- NULL

# read in RNA counts
# RNA <- read.csv("./results/data/SleuthGeneNormalized_TPM.csv", header=T) # Sleuth
RNA <- read.csv("./results/data/EdgeRGeneEstCountsMax.csv") # EdgeR
row.names(RNA) <- RNA$X
RNA$X <- NULL

# for sleuth
# colnames(RNA) <- c("TFAM_KO_Clone_11_1", "TFAM_KO_Clone_11_2", "TFAM_KO_Clone_5_1", "TFAM_KO_Clone_5_2", "TFAM_KO_Clone_6_1", "TFAM_KO_Clone_6_2", "HEK293T_NC.1_1", "HEK293T_NC.1_2", "HEK293T_NC.2_1", "HEK293T_NC.2_2", "HEK293T_NC.3_1", "HEK293T_NC.3_2")

# for edgeR
cols <- c("HEK293T_NC.1_1", "TFAM_KO_Clone_6_1", "HEK293T_NC.3_1", "TFAM_KO_Clone_11_1", "TFAM_KO_Clone_5_1", "HEK293T_NC.2_1", "HEK293T_NC.1_2", "TFAM_KO_Clone_6_2", "HEK293T_NC.3_2", "TFAM_KO_Clone_11_2", "TFAM_KO_Clone_5_2", "HEK293T_NC.2_2")
colnames(RNA) <- cols

RNA <- RNA[, colnames(EPIC)]

# Remove genes with TPM < 0.5 in >49% of samples
# RNA$TPM05 <- rowSums(RNA)
# RNA <- subset(RNA, RNA$TPM05 < 7)
# RNA$TPM05 <- NULL

# read in edgeR all genes to subset RNA counts by only genes from RNA analysis
genes <- read.csv("./results/data/EdgeR_RNA_all_genes.csv")
genes <- read.csv("./results/data/EdgeR_RNA_sig_genes.csv")

RNA <- subset(RNA, rownames(RNA) %in% genes$X)

# create dataframe to specify sample group Experiment (KO), Control (NC)
primary <- colnames(RNA)
GroupLabel <- c("Experiment", "Control", "Experiment", "Control", "Experiment", "Control", "Control", "Control", "Control", "Experiment", "Experiment", "Experiment")
Sample <- cbind(primary, GroupLabel)
Sample <- as.data.frame(Sample)
row.names(Sample) <- Sample$primary

# load EPIC manifest
sesameDataCache("EPIC.hg19.manifest")

# create MAE object for ELMER
data <- createMAE(exp = RNA, 
                  met = EPIC,
                  met.platform = "EPIC",
                  genome = "hg19",
                  save = FALSE,
                  TCGA = FALSE, 
                  colData=Sample,
                  met.na.cut=0.05
                  )

############################# ELMER find Probe-Gene Pairs ##########################

# find significant DMPs using ELMER
# not used in this analysis
# feeding signifcant from DMP analysis
# sig.diff <- get.diff.meth(data = data, 
#                           group.col = "GroupLabel",
#                           group1 =  "Control",
#                           group2 = "Experiment",
#                           mode = "supervised",
#                           minSubgroupFrac = 1, # if supervised mode set to 1
#                           sig.dif = 0,
#                           diff.dir = "both",
#                           cores = 1,  
#                           pvalue = 0.05,
#                           save = TRUE
#                           )

# read in DMP results to feed significant DMPs into ELMER
sig_cpgs <- as.data.frame(read.csv("./results/data/dmp_cont.csv"))
sig_cpgs <- filter(sig_cpgs, pval < 1e-7)

# find nearest 20 genes (10 upstream, 10 downstream) of significant CpGs
nearGenes <- GetNearGenes(data = data, 
                         probes = sig_cpgs$X, 
                         numFlankingGenes = 20)

# look at documentation: add diffExp = TRUE
pairs <- get.pair(data = data,
                      group.col = "GroupLabel",
                      group1 =  "Control",
                      group2 = "Experiment",
                      nearGenes = nearGenes,
                      mode = "supervised",
                      minSubgroupFrac = 1,
                      permu.size = 100, # Please set to 100000 to get significant results for unsupervised
                      raw.pvalue = 0.05,   
                      Pe = 0.05, # Please set to 0.001 to get significant results for unsupervised
                      diffExp = TRUE,
                      filter.probes = FALSE, # See preAssociationProbeFiltering function
                      filter.percentage = 0.05,
                      filter.portion = 0.3,
                      cores = 1,
                      label = "ALL_DMP_CONT_EDGER_MAX_5", 
                      diff.dir="both",
                      dir.out = './results/data')

save.image("./results/data/ELMER_TFAMKO_INTER_EDGER_MAX.RData")

pairs <- read.csv("./results/data/getPair.ALL_DMP_CONT_EDGER_MAX_SIG.pairs.statistic.with.empirical.pvalue.csv")
pairs1 <- read.csv("./results/data/getPair.ALL_DMP_CONT_EDGER_MAX_SIG.all.pairs.statistic.csv")
pairs2 <- read.csv("./results/data/getPair.ALL_DMP_CONT_EDGER_MAX_SIG.pairs.significant.csv")
pairs2 <- subset(pairs2, Experiment.vs.Control.diff.pvalue < 0.001)
pairs2 <- subset(pairs2, abs(Distance) < 1e6)
write.csv(pairs2,"./results/data/getPair.ALL_DMP_CONT_EDGER_MAX_SIG.pairs.significant.filtered.csv")
# log2FC_Experiment.vs.Control is the mean of the gene expression in the experiment vs mean of the gene expression in the control

enriched.motif <- get.enriched.motif(data = data,
                                     probes = pairs2$Probe, 
                                     label = "ALL_DMP_CONT_EDGER_MAX",
                                     min.incidence = 10,
                                     lower.OR = 1.1,
                                     save = TRUE,
                                     dir.out = './results/data')

TF <- get.TFs(data = data, 
              group.col = "GroupLabel",
              group1 =  "Experiment",
              group2 = "Control",
              mode = "supervised",
              enriched.motif = enriched.motif,
              cores = 1, 
              label = "ALL_DMP_CONT_EDGER_MAX", 
              diff.dir="both",
              dir.out = './results/data')

save.image("./results/data/ELMER_TFAMKO_FINAL_EDGER_MAX.RData")

################ ELMER Visualizations #################

# run the below block
CpG <- "cg05663891"
gene_id <- c("ENSG00000197134")
scatter.plot(data = data,
       byPair = list(probe = c(CpG), gene = gene_id), 
       category = "GroupLabel", save = TRUE, lm_line = TRUE, dir.out = "./results/plots/ELMER")

# need to again check if you see this in the other direction ... all 0's in TFAM and expression in controls

# Create scatterplots of each of Top 20 CpGs with 20 closest genes
probes20 <- sig_cpgs$X[1:20]
for (probes in probes20) {
  scatter.plot(data = data,
              byProbe = list(probe = c(probes), numFlankingGenes = 20), 
              category = "GroupLabel", 
              lm = TRUE, # Draw linear regression curve
              save = TRUE,
              dir.out = "./results/plots/ELMER") 
}

#schematic plot nearby genes
schematic.plot(pair = pairs2, 
               data = data,
               group.col = "GroupLabel",
               byProbe = pairs2$Probe[1],
               save = TRUE,
               dir.out = "./results/plots/ELMER")

#schematic plot nearby probes together
schematic.plot(pair = pairs, 
               data = data,   
               group.col = "GroupLabel", 
               byGene = pairs2$GeneID[1],
               save = TRUE,
               dir.out = "./results/plots/ELMER")

#TF expression vs. average DNA methylation try GABBR1
scatter.plot(data = data,
             byTF = list(TF = c("MAFK", "MAF", "MAFF", "MAFB", "MAFG"),
                         probe = enriched.motif[[names(enriched.motif)[1]]]), 
             category = "GroupLabel",
             save = TRUE, 
             lm_line = TRUE,
             dir.out = "./results/plots/ELMER")

# read in enrichment results to plot
enrich2 <- read.csv("getMotif.ALL_DMP_CONT_EDGER_MAX.motif.enrichment.csv", header=T)

motif.enrichment.plot(motif.enrichment = enrich2, 
                      significant = list(OR = 1.1,lowerOR = 1.1), 
                      label = "ALL_DMP_CONT_EDGER_MAX", 
                      summary = TRUE,
                      save = TRUE,
                      dir.out = "./results/plots/ELMER")  

heatmapPairs(data = data, 
             group.col = "GroupLabel",
             group1 =  "Experiment",
             group2 = "Control",
             pairs = pairs2,
             filename =  "HeatmapPairs.pdf")
             
# load("getTF.ALL.TFs.with.motif.pvalue.rda")
# motif <- colnames(TF.meth.cor)[1]
# TF.rank.plot(motif.pvalue = TF.meth.cor, 
#              motif = motif,
#              save = TRUE,
#              dir.out = ".results/plots/ELMER") 
        
# check how many genes in EdgeR overlap with ELMER
# 133 of 185 significant genes from RNA overlap with ELMER analysis
# elmer <- unique(pairs2$GeneID)
# edger <- read.csv("./results/data/EdgeR_RNA_sig_genes.csv")

# match <- logical(length(rem))
# for (i in 1:length(rem)){
#        for (j in 1:length(sig)){
#               if (rem[i] == sig[j]) {
#                      match[i] = TRUE
#               } 
#        }
# }
# table(match)