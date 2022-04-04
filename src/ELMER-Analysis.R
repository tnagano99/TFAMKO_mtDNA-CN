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
                      permu.size = 100000, # Please set to 100000 to get significant results
                      raw.pvalue = 0.05,   
                      Pe = 0.001, # Please set to 0.001 to get significant results
                     #  diffExp = TRUE,
                      filter.probes = FALSE, # See preAssociationProbeFiltering function
                      filter.percentage = 0.05,
                      filter.portion = 0.3,
                      cores = 1,
                      label = "ALL_DMP_CONT_EDGER_MAX", 
                      diff.dir="both",
                      dir.out = './results/data')

save.image("./results/data/ELMER_TFAMKO_INTER_EDGER_MAX.RData")

pairs <- read.csv("./results/data/getPair.ALL_DMP_CONT_EDGER_MAX.pairs.statistic.with.empirical.pvalue.csv")
pairs2 <- subset(pairs, pairs$FDR < 0.01)

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
CpG <- "cg05164933"
gene_id <- c("ENSG00000145423")
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

# match <- logical(length(edger$X))
# for (i in 1:length(edger$X)){
#        for (j in 1:length(elmer)){
#               if (edger$X[i] == elmer[j]) {
#                      match[i] = TRUE
#               } 
#        }
# }
# table(match)

get.pair <- function(data,
                     nearGenes,
                     minSubgroupFrac = 0.4,
                     permu.size = 10000,
                     permu.dir = NULL,
                     raw.pvalue = 0.001,
                     Pe = 0.001,
                     mode = "unsupervised",
                     diff.dir = NULL,
                     dir.out = "./",
                     diffExp = FALSE,
                     group.col,
                     group1 = NULL,
                     group2 = NULL,
                     cores = 1,
                     correlation = "negative",
                     filter.probes = TRUE,
                     filter.portion = 0.3,
                     filter.percentage = 0.05,
                     label = NULL,
                     addDistNearestTSS = FALSE,
                     save = TRUE){
  
  if(is.character(nearGenes)){
    nearGenes <- get(load(nearGenes))
  }
  
  # if(!all(c("ID", "GeneID", "Symbol" ) %in% colnames(nearGenes)))
  #   stop("nearGenes does not have one of the expected columns: ID, GeneID, Symbol")
  
  if(diffExp & missing(group.col))
    stop("Please set group.col argument to test whether putative target gene are differentially expressed between two groups.")
  
  if(missing(group.col)) stop("Please set group.col argument")
  if(missing(group1)) stop("Please set group1 argument")
  if(missing(group2)) stop("Please set group2 argument")
  data <- data[,colData(data)[,group.col] %in% c(group1, group2)]
  
  # Supervised groups
  unmethylated <- methylated <- NULL
  if(mode == "supervised"){
    if(is.null(diff.dir)) stop("For supervised mode please set diff.dir argument (same from the get.diff.meth)")
    if(diff.dir == "hypo"){
      message("Using pre-defined groups. U (unmethylated): ",group1,", M (methylated): ", group2)
      unmethylated <-  which(colData(data)[,group.col]  == group1)
      methylated <-  which(colData(data)[,group.col]  == group2)
    } else {
      message("Using pre-defined groups. U (unmethylated): ",group2,", M (methylated): ", group1)
      unmethylated <-  which(colData(data)[,group.col]  == group2)
      methylated <-  which(colData(data)[,group.col]  == group1)
    }
  } else {
    message("Selecting U (unmethylated) and M (methylated) groups. Each groups has ", minSubgroupFrac * 50,"% of samples")
  }
  # Paralellization code
  parallel <- FALSE
  if (cores > 1){
    if (cores > detectCores()) cores <- detectCores()
    registerDoParallel(cores)
    parallel = TRUE
  }
  
  if(filter.probes) data <- preAssociationProbeFiltering(data, K = filter.portion, percentage = filter.percentage)
  
  met <- assay(getMet(data))
  # Probes that were removed from the last steps cannot be verified
  nearGenes <- nearGenes[nearGenes$ID %in% rownames(met),]
  
  if(nrow(nearGenes) == 0) {
    message("No probes passed the preAssociationProbeFiltering filter")
    return(NULL)
  }
  exp <- assay(getExp(data))
  message("Calculating Pp (probe - gene) for all nearby genes")
  Probe.gene <- adply(.data = unique(nearGenes$ID),
                      .margins = 1,
                      .fun = function(x) {
                        Stat.nonpara(
                          Probe = x,
                          Meths = met[x,],
                          methy = methylated,
                          unmethy = unmethylated,
                          NearGenes = as.data.frame(nearGenes),
                          correlation = correlation,
                          Top = minSubgroupFrac/2, # Each group will have half of the samples
                          Exps = exp
                        )
                      },
                      .progress = "time",
                      .parallel = parallel,
                      .id = NULL,
                      .paropts = list(.errorhandling = 'pass')
  )
  
  rownames(Probe.gene) <- paste0(Probe.gene$Probe,".",Probe.gene$GeneID)
  Probe.gene <- Probe.gene[!is.na(Probe.gene$Raw.p),]
  
  if(save) {
    dir.create(dir.out, showWarnings = FALSE)
    file <- sprintf("%s/getPair.%s.all.pairs.statistic.csv",dir.out, ifelse(is.null(label),"",label))
    write_csv(Probe.gene,file = file)
    message(paste("File created:", file))
  }
  
  Probe.gene <- Probe.gene[Probe.gene$Raw.p < raw.pvalue,]
  Probe.gene <- Probe.gene[order(Probe.gene$Raw.p),]
  selected <- Probe.gene
  if(nrow(selected) == 0) {
    message(paste("No significant pairs were found for pvalue =", raw.pvalue))
    return(selected)
  }
  
  
  #   Probe.gene$logRaw.p <- -log10(Probe.gene$Raw.p)
  if(mode == "unsupervised"){
    GeneID <- unique(Probe.gene[,"GeneID"])
    message(paste("Calculating Pr (random probe - gene). Permutating ", permu.size, "probes for",  length(GeneID), "genes"))
    # get permutation
    permu <- get.permu(data,
                       geneID     = GeneID,
                       percentage = minSubgroupFrac / 2,
                       rm.probes  = unique(nearGenes$ID),
                       methy      = methylated,
                       unmethy    = unmethylated,
                       correlation = correlation,
                       permu.size = permu.size,
                       permu.dir  = permu.dir,
                       cores      = cores)
    # Get empirical p-value
    Probe.gene.Pe <- Get.Pvalue.p(Probe.gene,permu)
    
    if(save) write_csv(Probe.gene.Pe,
                       file = sprintf("%s/getPair.%s.pairs.statistic.with.empirical.pvalue.csv",dir.out,
                                      ifelse(is.null(label),"",label)))
    # Pe will always be 1 for the supervised mode. As the test exp(U) > exp(M) will always be doing the same comparison.
    selected <- Probe.gene.Pe[Probe.gene.Pe$Pe < Pe & !is.na(Probe.gene.Pe$Pe),]
  } else {
    Probe.gene$FDR <- p.adjust(Probe.gene$Raw.p,method = "BH")
    if(save) write_csv(Probe.gene,
                       file=sprintf("%s/getPair.%s.pairs.statistic.with.empirical.pvalue.csv",dir.out,
                                    ifelse(is.null(label),"",label)))
    selected <-  Probe.gene[Probe.gene$FDR < Pe,]
  }
  
  # Change distance from gene to nearest TSS
  if(addDistNearestTSS) {
    selected$Distance <- NULL
    selected <- addDistNearestTSS(data, NearGenes = selected)
  }
  if(diffExp){
    message("Calculating differential expression between two groups")
    Exp <- assay(getExp(data)[unique(selected$GeneID),])
    groups <- unique(colData(data)[,group.col])
    prefix <- paste(gsub("[[:punct:]]| ", ".", groups),collapse =  ".vs.")
    log.col <- paste0("log2FC_",prefix)
    diff.col <- paste0(prefix,".diff.pvalue")
    idx1 <- colData(data)[,group.col] == groups[1]
    idx2 <- colData(data)[,group.col] == groups[2]
    out <- adply(.data = split(Exp,rownames(Exp)), .margins = 1,
                 .fun = function(x) {
                   test <- t.test(x = x[idx1],y = x[idx2])
                   out <- data.frame("log2FC" = test$estimate[1] - test$estimate[2],
                                     "diff.pvalue" = test$p.value)
                 },
                 .progress = "time",
                 .parallel = parallel,
                 .id = "GeneID",
                 .paropts = list(.errorhandling = 'pass')
    )
    head(add)
    add <- out[match(selected$GeneID, out$GeneID),c("log2FC","diff.pvalue")]
    colnames(add) <- c(log.col,diff.col)
    selected <- cbind(selected, add)
  }
  if(save) write_csv(selected,file = sprintf("%s/getPair.%s.pairs.significant.csv",dir.out, ifelse(is.null(label),"",label)))
  invisible(gc())
  return(selected)
}