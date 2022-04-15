#' chromHMM
#' @title Chromatin state enrichments
#' @description Perfom enrichment testing against the ChromHMM 15-state model for hg38
#'  using \code{LOLA}
#' @param sigRegions A \code{GRanges} object of significant regions
#' @param regions A \code{GRanges} object of background regions
#' @param cores An integer of how many cores to use
#' @return A \code{tibble} of enrichment results
#' @importFrom dplyr as_tibble select mutate summarize pull mutate_if arrange recode_factor
#' @importFrom tidyr pivot_wider
#' @importFrom LOLA loadRegionDB runLOLA writeCombinedEnrichment
#' @importFrom magrittr %>% %T>%
#' @importFrom hablar s
#' @export chromHMM
#' 
chromHMM <- function(sigRegions = sigRegions,
                     regions = regions,
                     cores = cores){
  message("Performing ChromHMM enrichment testing")
  chromHMM <- LOLA::loadRegionDB(dbLocation = "./DMRichR/all.mnemonics.bedFiles",
                                 useCache = TRUE,
                                 limit = NULL,
                                 collections = "Roadmap_ChromHMM") %>%
    LOLA::runLOLA(userSets = sigRegions,
                  userUniverse = regions,
                  regionDB = .,
                  minOverlap = 1,
                  cores = cores,
                  redefineUserSets = FALSE) %T>%
    LOLA::writeCombinedEnrichment(combinedResults = .,
                                  outFolder = "ChromHMM",
                                  includeSplits = FALSE) %>%
    dplyr::as_tibble() %>%
    dplyr::select(oddsRatio, cellType, tissue, antibody) %>%
    dplyr::mutate(antibody = as.factor(antibody)) %>% 
    dplyr::mutate(antibody = dplyr::recode_factor(antibody,
                                                  "1_TssA" = "Active TSS",
                                                  "2_TssAFlnk" = "Flanking Active TSS",
                                                  "3_TxFlnk" = "Transcription at Gene 5' and 3'",
                                                  "4_Tx" = "Strong Transcription",
                                                  "5_TxWk" = "Weak Transcription",
                                                  "6_EnhG"= "Genic Enhancers",
                                                  "7_Enh" = "Enhancers",
                                                  "8_ZNF/Rpts" = "ZNF Genes & Repeats",
                                                  "9_Het" = "Heterochromatin",
                                                  "10_TssBiv" = "Bivalent/Poised TSS",
                                                  "11_BivFlnk" = "Flanking Bivalent TSS/Enhancer",
                                                  "12_EnhBiv" = "Bivalent Enhancer",
                                                  "13_ReprPC" = "Repressed PolyComb",
                                                  "14_ReprPCWk" = "Weak Repressed PolyComb",
                                                  "15_Quies" = "Quiescent/Low"
    )
    ) %>%
    dplyr::arrange(antibody) 
  
  # Fix Inf Odds Ratio
  max <- chromHMM %>%
    dplyr::summarize(max(hablar::s(oddsRatio))) %>% 
    dplyr::pull()
  
  chromHMM <- chromHMM %>%
    dplyr::mutate_if(is.numeric, function(x) ifelse(is.infinite(x), max, x)) %>% 
    tidyr::pivot_wider(names_from = antibody, values_from = oddsRatio) %>%
    dplyr::arrange(tissue) %>%
    return()
}

#' chromHMM_heatmap
#' @title Chromatin state heatmap
#' @description Plot a heatmap of \code{LOLA} enrichment testing results of the ChromHMM
#'  15-state model for hg38
#' @param chromHMM A \code{tibble} of enrichment results
#' @return Saves a heatmap
#' @importFrom dplyr group_by tally select
#' @importFrom magrittr %>%
#' @importFrom viridis viridis
#' @importFrom gplots heatmap.2
#' @importFrom PerformanceAnalytics tol21rainbow
#' @export chromHMM_heatmap
#' 
chromHMM_heatmap <- function(chromHMM = chromHMM){
  message("Plotting ChromHMM heatmap")
  
  # Make Row Labels
  labels <- chromHMM %>%
    dplyr::group_by(tissue) %>% 
    dplyr::tally()
  
  # Colors
  palette(PerformanceAnalytics::tol21rainbow)
  rowcolors <- as.list(palette(PerformanceAnalytics::tol21rainbow))
  rowcolors <- rowcolors[1:nrow(labels)]
  
  colorlist <- list()
  for (i in 1:length(labels$n)){
    colorlist[[i]] <- rep(rowcolors[i], labels$n[i])
  }
  colorlist <- unlist(colorlist)
  
  # Select matrix data
  data <- chromHMM %>%
    dplyr::select(-cellType, -tissue) %>%
    data.matrix()
  
  # Plot Heatmap
  pdf("ChromHMM/ChromHMM_heatmap.pdf",
      height = 8.5,
      width = 12)
  
  gplots::heatmap.2(data,
                    Rowv = F,
                    Colv = F,
                    dendrogram = "none",
                    col = viridis::viridis(15, option = "inferno"),
                    margins = c(15,2),
                    trace = "none",
                    labRow = "" ,
                    labCol = colnames(data),
                    main = "Enriched Chromatin States for Differentially Methylated Regions",
                    RowSideColors = colorlist,
                    srtCol = 60,
                    keysize = 0.85,
                    key.par = list(cex = 0.5),
                    key.xlab = "Odds Ratio",
                    key.ylab = "Frequency",
                    key.title = ""
  )
  
  # Legend
  par(xpd = TRUE, mar = par()$mar + c(0,6,0,0))
  legend(x = -0.075,
         y= 0.9,
         legend = labels$tissue,
         col = unlist(rowcolors),
         lty = 1,
         lwd = 6,
         cex = 1,
         bty = "n")
  dev.off()
}

#' roadmap
#' @title Chromatin mark enrichments
#' @description Perfom enrichment testing against the Roadmap Epigenomics core histone modifications
#'  (5 marks, 127 epigenomes) for hg38 using \code{LOLA}
#' @param sigRegions A \code{GRanges} object of significant regions
#' @param regions A \code{GRanges} object of background regions
#' @param cores An integer of how many cores to use
#' @return A \code{tibble} of enrichment results
#' @importFrom dplyr as_tibble select mutate rename summarize pull mutate_if arrange left_join 
#'  group_by tally filter
#' @importFrom tidyr pivot_wider replace_na
#' @importFrom stringr str_replace str_to_title
#' @importFrom LOLA loadRegionDB runLOLA writeCombinedEnrichment
#' @importFrom magrittr %>% %T>%
#' @importFrom hablar s
#' @export roadmap
#' 
roadmap <- function(sigRegions = sigRegions,
                    regions = regions,
                    cores = cores){
  message("Performing Roadmap epigenomics enrichment testing")
  # update dbLocation to hg19 folder in data/roadmap_epigenomics/hg19
  roadmap <- LOLA::loadRegionDB(dbLocation = "./DMRichR/LOLARoadmap_180423/hg19",
                                useCache = TRUE, 
                                limit = NULL,
                                collections = "roadmap_epigenomics") %>%
    LOLA::runLOLA(userSets = sigRegions,
                  userUniverse = regions,
                  regionDB = .,
                  minOverlap = 1,
                  cores = cores,
                  redefineUserSets = FALSE) %T>%
    LOLA::writeCombinedEnrichment(combinedResults = .,
                                  outFolder = "RoadmapEpigenomics",
                                  includeSplits = FALSE) %>%
    dplyr::as_tibble() %>%
    # update hg19 folder in data/roadmap_epigenomics/hg19
    dplyr::left_join(readr::read_tsv("./DMRichR/LOLARoadmap_180423/hg19/roadmap_epigenomics/index.txt"),
                     by = "filename") %>%
    dplyr::select(oddsRatio, antibody.x, donor, anatomy) %>%
    dplyr::rename(antibody = antibody.x) %>% 
    dplyr::mutate(antibody = tidyr::replace_na(antibody, "DNase"))
  
  # Subset for Core Histone Mods from same 127 samples: H3K4me1, H3K4me3, H3K27me3, H3K36me3, H3K9me3
  core <- roadmap %>%
    dplyr::group_by(antibody) %>% 
    dplyr::tally() %>%
    dplyr::filter(n == 127) %>%
    dplyr::pull(antibody)
  
  roadmap <- roadmap %>% 
    dplyr::filter(antibody %in% core) %>%
    dplyr::arrange(anatomy)
  
  # Fix Inf Odds Ratio
  max <- roadmap %>%
    dplyr::summarize(max(hablar::s(oddsRatio))) %>% 
    dplyr::pull()
  
  roadmap <- roadmap %>%
    dplyr::mutate_if(is.numeric, function(x) ifelse(is.infinite(x), max, x)) %>% 
    tidyr::pivot_wider(names_from = antibody, values_from = oddsRatio) %>%
    dplyr::mutate(anatomy = stringr::str_replace(anatomy, "_", " ")) %>%
    dplyr::mutate(anatomy = stringr::str_to_title(anatomy)) %>%
    dplyr::mutate(anatomy = stringr::str_replace(anatomy, "Ipsc", "IPSC")) %>%
    dplyr::mutate(anatomy = stringr::str_replace(anatomy, "Esc", "ESC")) %>%
    dplyr::mutate(anatomy = stringr::str_replace(anatomy, "Gi", "GI")) %>%
    return()
}

#' roadmap_heatmap
#' @title Chromatin mark heatmap
#' @description Plot a heatmap of \code{LOLA} enrichment testing results of the Roadmap Epigenomics
#'  core marks for hg38
#' @param roadmap A \code{tibble} of enrichment results
#' @return Saves a heatmap
#' @importFrom dplyr group_by tally select
#' @importFrom magrittr %>%
#' @importFrom viridis viridis
#' @importFrom gplots heatmap.2
#' @importFrom PerformanceAnalytics tol21rainbow
#' @export roadmap_heatmap
roadmap_heatmap <- function(roadmap = roadmap){
  message("Plotting Roadmap Epigenomics heatmap")
  
  # Make Row Labels
  labels <- roadmap %>%
    dplyr::group_by(anatomy) %>% 
    dplyr::tally() 
  
  # Colors
  palette(PerformanceAnalytics::tol21rainbow)
  rowcolors <- as.list(palette(PerformanceAnalytics::tol21rainbow))
  rowcolors <- rep(rowcolors, 2)
  rowcolors <-rowcolors[1:nrow(labels)]
  
  colorlist <- list()
  for (i in 1:length(labels$n)){
    colorlist[[i]] <- rep(rowcolors[i], labels$n[i])
  }
  colorlist <- unlist(colorlist)
  
  # Select matrix data
  data <- roadmap %>%
    dplyr::select(-donor, -anatomy) %>%
    data.matrix()
  
  # Plot Heatmap
  pdf("RoadmapEpigenomics/Roadmap_heatmap.pdf",
      height = 8.5,
      width = 12)
  
  gplots::heatmap.2(data,
                    Rowv = F,
                    Colv = F,
                    dendrogram = "none",
                    # col = viridis::viridis(15, option = "inferno"),
                    col = rev(heat.colors(15, 1)),
                    margins =c (15,2),
                    trace = "none",
                    labRow = "" ,
                    labCol = colnames(data),
                    main = "Histone Modifications at Differentially Methylated Regions",
                    RowSideColors = colorlist,
                    srtCol = 60,
                    keysize = 0.85,
                    key.par = list(cex = 0.5),
                    key.xlab = "Odds Ratio",
                    key.ylab = "Frequency",
                    key.title = "",
                    density.info = 'none'
  )
  
  # Legend
  par(xpd = TRUE, mar = par()$mar + c(0,6,0,0))
  legend(x = -0.075,
         y= 0.9,
         legend = unlist(labels$anatomy),
         col = unlist(rowcolors),
         lty = 1,
         lwd = 6,
         cex = 1,
         bty = "n")
  dev.off()
}

#' dmrList
#' @title Stratify DMRs by directionality 
#' @description Create \code{GRangesList} object of all, hypermethylated,
#'  and hypomethylated DMRs from \code{dmrseq::dmrseq()}
#' @param sigRegions A \code{GRanges} object of DMRs from \code{dmrseq::dmrseq()}
#' @return A \code{GRangesList} of DMRs
#' @importFrom magrittr %>%
#' @importFrom plyranges filter
#' @importFrom GenomicRanges GRangesList
#' @export dmrList
#' 
dmrList <- function(sigRegions = sigRegions){
  message("Making DMR list")
  
  GenomicRanges::GRangesList("All DMRs" = sigRegions,
                             "Hypermethylated DMRs" = sigRegions %>%
                               plyranges::filter(stat > 0),
                             "Hypomethylated DMRs" = sigRegions %>%
                               plyranges::filter(stat < 0)) %>% 
    return()
}

##################################################################################

########################### Roadmap Epigenomics #################################
library(dplyr)
library(tidyr)
library(LOLA)
library(hablar)
library(simpleCache)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(data.table)
library(GenomicRanges)
library(magrittr)
library(glue)
library(stringr)
library(GenomeInfoDb)
library(forcats)
library(plyranges)
library(stats)
library(DMRichR)
library(ggplot2)

ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)


## Load in DMPs --------------------------------------------------------------
# update to dmp_cont in results/data
dmps <- read.csv("./DMRichR/dmp_cont.csv")
tmp <- dmps$X
rownames(dmps) <- dmps$X
dmps$X <- NULL
dmps[] <- sapply(dmps, as.double)
dmps$X <- tmp

## Create bed files ------------------------------------------------------------
StatsAnn <- merge(dmps, ann, by.x="X", by.y="row.names")
StatsAnn$pos2 <- StatsAnn$pos

# REGION UNIVERSE: Bed file of chromosome start, stop, name
statsbed <- StatsAnn[, c("chr", "pos", "pos2", "X")]
statsbed <- as.data.frame(statsbed)

## QUERY SET: lists of genomic regions to be tested for enrichment
sigStats <- statsbed[StatsAnn$pval < 1e-7,] 


## Turning dataframe into genomic ranges ---------------------------------------
statsbed <- makeGRangesFromDataFrame(statsbed, 
                                     start.field="pos", end.field="pos2", ignore.strand=T)
sigStats <- makeGRangesFromDataFrame(sigStats, 
                                     start.field="pos", end.field="pos2", ignore.strand=T)

roadmapEnrichments <- roadmap(sigStats, statsbed, cores = 1)

####################### ChromHMM ###############################################

chromHMMEnrichments <- chromHMM(sigStats, statsbed, cores = 1)

####################### Data prep for chromHMM #################################

# unzip all gz files to get bed files
# set working directory to all gz files
library(R.utils)

filenames <- dir()
for (i in 1:length(filenames)){
  gunzip(filenames[i])
}

# read in each bed file and append to dataframe
filenames <- dir()
df <- data.frame(matrix(ncol = 4, nrow = 0))

for (i in 1:length(filenames)){
  epi <- as.data.frame(read.table(filenames[i], header = FALSE, sep = "\t", quote=""))
  df <- rbind(df, epi)
}
colnames(df) <- c("chr", "start", "end", "group")

# create subsets of each group

groups <- c("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF/Rpts","9_Het",
            "10_TssBiv","11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","15_Quies")

ret <- GRangesList()
size <- numeric(length = length(groups))
for (i in 1:length(groups)){
  cat(i)
  epi <- df[df$group == groups[i],]
  size[i] <- nrow(epi)
  epibed <- makeGRangesFromDataFrame(epi, start.field="start", end.field="end")
  ret <- c(ret, GRangesList(epibed))
}
save(ret, file = "Roadmap_ChromHMM.Rdata")

# use data to create annotation sheet
df_anno <- data.frame(matrix(ncol = 9, nrow = 15))
colnames(df_anno) <- c("filename","cellType","description","tissue","dataSource","antibody","treatment","collection","size")

df_anno$filename <- groups
df_anno$collection <- "Roadmap_ChromHMM"
df_anno$size <- size
df_anno$antibody <- groups
df_anno$description <- paste(df_anno$collection, df_anno$antibody, sep = " ")
df_anno$cellType <- as.character(df_anno$cellType)
df_anno$tissue <- as.character(df_anno$tissue)
df_anno$dataSource <- as.character(df_anno$dataSource)
df_anno$treatment <- as.character(df_anno$treatment)

ret <- df_anno
save(ret, file = "Roadmap_ChromHMM_files.RData")
