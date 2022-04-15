## DMRichR Functions by Ben Laufer ---------------------------------------------
#  https://github.com/ben-laufer/DMRichR

#' DMRichGenic
#' @title Gene region enrichment testing
#' @description Test DMRs for overlaps with gene region annotations for all genomes
#' @param sigRegions \code{GRanges} object of DMRs
#' @param regions \code{GRanges} object of background regions
#' @param annoDb Character specifying OrgDb annotation package for species of interest
#' @param TxDb TxDb annotation package for genome of interest
#' @return A tibble with the enrichment results
#' @importFrom dplyr filter mutate case_when select as_tibble
#' @importFrom magrittr %>%
#' @importFrom glue glue glue_collapse
#' @importFrom stringr str_detect
#' @importFrom data.table rbindlist
#' @importFrom GenomeInfoDb genome
#' @importFrom forcats as_factor
#' @importFrom plyranges as_granges
#' @importFrom stats fisher.test p.adjust
#' @export DMRichGenic

DMRichGenic <- function(sigRegions = sigRegions,
                        regions = regions,
                        TxDb = TxDb,
                        annoDb = annoDb){
  
  genome <- TxDb %>%
    GenomeInfoDb::genome() %>%
    unique()
  
  sigRegions <- sigRegions %>%
    plyranges::as_granges()
  
  regions <- regions %>% 
    plyranges::as_granges()
  
  print(glue::glue("{genome} annotations will be used for sigRegions and regions"))
  sigRegionsAnnotated <- sigRegions %>%
    DMRichR::annotateRegions(TxDb,
                             annoDb)
  
  regionsAnnotated <- regions %>%
    DMRichR::annotateRegions(TxDb,
                             annoDb)
  
  annotations <- c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Downstream", "Intergenic")
  print(glue::glue("Performing gene region enrichment testing for {genome} with the following annotations: {tidyAnnotations}",
                   tidyAnnotations = glue::glue_collapse({annotations}, sep = ", ", last = ", and ")))
  
  lapply(annotations, function(term){
    
    print(glue::glue("Now FISHERing for {term} annotation"))
    
    sigRegionsOverlap <- sigRegionsAnnotated %>%
      dplyr::filter(stringr::str_detect(annotation, term))
    
    regionsOverlap <- regionsAnnotated %>%
      dplyr::filter(stringr::str_detect(annotation, term))
    
    genicMatrix <- matrix(c(nrow(sigRegionsOverlap), (length(sigRegions) - nrow(sigRegionsOverlap)),
                            nrow(regionsOverlap), (length(regions) - nrow(regionsOverlap))),
                          nrow = 2)
    
    results <- data.frame("Annotation" =  term,
                          "OR" = fisher.test(genicMatrix)[["estimate"]][["odds ratio"]],
                          "CIlower" = fisher.test(genicMatrix)[["conf.int"]][1],
                          "CIupper" = fisher.test(genicMatrix)[["conf.int"]][2],
                          "p.value" = fisher.test(genicMatrix)$p.value)
    
  }) %>%
    data.table::rbindlist() %>%
    dplyr::mutate(fdr = p.adjust(p.value, method = 'fdr')) %>%
    dplyr::mutate(Annotation = forcats::as_factor(Annotation)) %>% 
    dplyr::as_tibble() %>% 
    return()
}

#' DMRichCpG
#' @title CpG annotation enrichment testing
#' @description Test DMRs for overlaps with CpG annotations (islands, shores, shelves, open sea)
#' @param sigRegions \code{GRanges} object of DMRs
#' @param regions \code{GRanges} object of background regions
#' @param genome A character vector specifying the genome of interest
#'  c("hg38", "hg19", "mm10", "mm9", "rheMac10", "rheMac8", "rn6", "danRer11",
#'   "galGal6", "bosTau9", "panTro6", "dm6", "susScr11", "canFam3")
#' @return A tibble with the enrichment results
#' @importFrom dplyr filter mutate case_when select recode_factor as_tibble
#' @importFrom magrittr %>%
#' @importFrom glue glue glue_collapse
#' @importFrom stringr str_detect
#' @importFrom data.table rbindlist
#' @importFrom GenomeInfoDb keepStandardChromosomes genome
#' @importFrom plyranges as_granges mutate count_overlaps
#' @importFrom stats fisher.test p.adjust
#' @export DMRichCpG

DMRichCpG <- function(sigRegions = sigRegions,
                      regions = regions,
                      genome = genome){
  
  stopifnot(genome %in% c("hg38", "hg19", "mm10", "mm9", "rheMac10", "rheMac8", "rn6", "danRer11", "galGal6", "bosTau9", "panTro6", "dm6", "susScr11", "canFam3"))
  print(glue::glue("Performing CpG annotation enrichment testing for {genome}"))
  
  CGannotations <- genome %>%
    DMRichR::getCpGs()
  
  sigRegions <- sigRegions %>%
    plyranges::as_granges()
  
  regions <- regions %>% 
    plyranges::as_granges()
  
  CGterms <- c("islands", "shores", "shelves", "inter")
  
  lapply(CGterms, function(term){
    
    print(glue::glue("Now FISHERing for {term} annotation"))
    
    CGannotationsFiltered <- CGannotations %>%
      plyranges::filter(stringr::str_detect(type, term))
    
    sigRegionsOverlap <- sigRegions %>% 
      plyranges::mutate(n_overlaps = plyranges::count_overlaps(.,CGannotationsFiltered)) %>%
      plyranges::filter(n_overlaps > 0)
    
    regionsOverlap <- regions %>% 
      plyranges::mutate(n_overlaps = plyranges::count_overlaps(.,CGannotationsFiltered)) %>%
      plyranges::filter(n_overlaps > 0)
    
    CGmatrix <- matrix(c(length(sigRegionsOverlap), (length(sigRegions) - length(sigRegionsOverlap)),
                         length(regionsOverlap), (length(regions) - length(regionsOverlap))),
                       nrow = 2)
    
    results <- data.frame("Annotation" =  term,
                          "OR" = fisher.test(CGmatrix)[["estimate"]][["odds ratio"]],
                          "CIlower" = fisher.test(CGmatrix)[["conf.int"]][1],
                          "CIupper" = fisher.test(CGmatrix)[["conf.int"]][2],
                          "p.value" = fisher.test(CGmatrix)$p.value)
    
  }) %>%
    data.table::rbindlist() %>%
    dplyr::mutate(fdr = p.adjust(p.value, method = 'fdr')) %>%
    dplyr::mutate(Annotation = dplyr::recode_factor(Annotation,
                                                    islands = "CpG Islands",
                                                    shores = "CpG Shores",
                                                    shelves = "CpG Shelves",
                                                    inter = "Open Sea")
    ) %>% 
    dplyr::as_tibble() %>% 
    return()
}

#' DMRichPlot
#' @title Plot gene region and CpG annotation enrichment testing
#' @description Plot DMR CpG or gene region enrichment testing results
#'  from \code{DMRichR::DMRichGenic()} or \code{DMRichR::DMRichCpG()}.
#' @param data A \code{tibble} from \code{DMRichR::DMRichGenic()} or \code{DMRichR::DMRichCpG()}.
#' @param type A character vector of the type of results to plot i.e. c("CpG", "genic").
#' @param multi A logical indicating whether to create facets from \code{DMRichR::DMparseR)} output.
#' @return A \code{ggplot} object of enrichment results that can be viewed by calling it, 
#' saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @import ggplot2
#' @importFrom dplyr select mutate
#' @importFrom forcats as_factor
#' @importFrom wesanderson wes_palette
#' @export DMRichPlot
#' 

DMRichPlot <- function(data = data,
                       type = c("CpG", "genic"),
                       multi = FALSE){
  
  stopifnot(type %in% c("CpG", "genic"))
  print(glue::glue("Plotting {type} annotation results"), "\n")
  
  ########## added between
  
  # updates x-axis label depending on CpG or Genic
  if (type == "CpG") {
    x_label <- "CpG Functional Groups"
  } else {
    x_label <- "Functional Genomic Elements"
  }
  
  # sets plot title
  plot_title <- "Functional Region Enrichments"
  
  ##########
  
  data <- data %>%
    dplyr::mutate(OR = dplyr::case_when(OR < 1 ~ -1/OR,
                                        OR >= 1 ~ OR)
    ) %>%
    dplyr::mutate(signif = dplyr::case_when(fdr <= 0.05 ~ 1,
                                            fdr > 0.05 ~ 0)
    ) 
  
  p <- ggplot(data = data,
              aes(x = Annotation,
                  y = OR,
                  fill = -log(fdr, base = 10)) # updated fill from Annotation
  ) +
    geom_bar(stat = "identity", 
             color = "Black") +
    coord_flip() +
    labs(y = "Fold Enrichment",
         x = x_label, # updated depending on CpG or genic
         title = plot_title, # updated depending on CpG or genic
         fill = "-log10(FDR)"
    ) +
    theme_classic() + 
    theme(plot.title = element_text(size = 24, hjust = 0.5, face = "bold"), # added to set plot title size and center
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.position = "right" # updated from "none"
    ) +
    scale_y_continuous(expand = c(0.1, 0.1)) + 
    scale_x_discrete(limits = data$Annotation %>%
                       as.factor() %>% 
                       levels() %>%
                       rev()
    ) + 
    geom_hline(yintercept = 0) +
    geom_text(data = data[(data$signif == 1 & data$OR > 0), ],
              label = "*",
              size = 8,
              show.legend = FALSE,
              nudge_y = 0.1, # updated to move star closer to bar
              nudge_x = -0.09) +
    geom_text(data = data[(data$signif == 1 & data$OR < 0), ],
              label = "*",
              size = 8,
              show.legend = FALSE,
              nudge_y = -0.1, # updated to move star closer to bar
              nudge_x = -0.09)
  
  if(type == "CpG"){
    p <- p +
      scale_fill_manual(values = c("forestgreen", "goldenrod2", "dodgerblue", "blue3"),
                        breaks = data$Annotation %>%
                          forcats::as_factor() %>% 
                          levels(),
                        name = "Annotation")
  }else if(type == "genic"){
    # updated this block to use scale_fill_gradient
    # p <- p +
    #   scale_fill_manual(values = data$Annotation %>%
    #                       as.factor() %>%
    #                       nlevels() %>%
    #                       wesanderson::wes_palette("Zissou1", n = ., type = "continuous") %>%
    #                       rev(),
    #                     breaks = data$Annotation %>%
    #                       levels(),
    #                     name = "Annotation")
    
    p <- p +
      scale_fill_gradient(low='white', high='red', space='Lab')
  }
  if(multi == TRUE){
    p <- p +
      facet_grid(~Direction)
  }
  return(p)
}

#' DMparseR
#' @title Combined plot of stratified gene region and CpG annotation enrichment testing
#' @description Parse the results from multiple \code{DMRichR} runs.
#' @param direction A character vector of the DMR profiles to analyze
#'  c("All DMRs", "Hypermethylated DMRs", "Hypomethylated DMRs").
#' @param type A character vector of the type of results to parse i.e. c("CpG", "genic").
#' @return A \code{tibble} of enrichment results
#' @importFrom dplyr filter mutate case_when select as_tibble
#' @importFrom tidyr separate
#' @importFrom magrittr %>%
#' @importFrom glue glue glue_collapse
#' @importFrom data.table rbindlist
#' @importFrom readxl read_xlsx
#' @importFrom forcats as_factor
#' @export DMparseR
#' 
DMparseR <- function(direction = c("All DMRs", "Hypermethylated DMRs", "Hypomethylated DMRs"),
                     type = c("CpG", "genic")
){
  
  stopifnot(direction %in% c("All DMRs", "Hypermethylated DMRs", "Hypomethylated DMRs"))
  stopifnot(type %in% c("CpG", "genic"))
  
  print(glue::glue("Parsing {type} enrichment results for {tidyDirection}",
                   tidyDirection = glue::glue_collapse({direction}, sep = ", ", last = " and ")))
  
  purrr::map(direction,
             function(direction){
               glue::glue("DMRichments/{direction}_{type}_enrichments.xlsx") 
             }) %>%
    as.vector() %>%
    lapply(function(file){
      readxl::read_xlsx(file)
    }) %>%
    `names<-` (direction) %>% 
    data.table::rbindlist(idcol = "Dataset") %>%
    dplyr::as_tibble() %>%
    tidyr::separate(Dataset, c("Direction", "DMR")) %>%
    dplyr::select(Direction,
                  Annotation,
                  OR,
                  fdr) %>%
    dplyr::mutate(Annotation = forcats::as_factor(Annotation)) %>%
    return()
}





## Load in packages ------------------------------------------------------------
library(simpleCache)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(data.table)
library(GenomicRanges)
library(dplyr)
library(magrittr)
library(glue)
library(stringr)
library(GenomeInfoDb)
library(forcats)
library(plyranges)
library(stats)
library(DMRichR)
library(ggplot2)
library(wesanderson)
library(xlsx)

# Load in annotation file
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)


## Load in DMPs --------------------------------------------------------------
setwd('C:/Users/tnaga/Documents/R_Scripts/Thesis/DMRichR')

dmps <- read.csv("dmp_cont.csv")
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


## Stratified Gene Region Enrichment Test --------------------------------------
# significant results have and fdr < 0.05
RegionEnrichment <- DMRichGenic(sigStats, statsbed, 
                                TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                "org.Hs.eg.db")

Region_plot <- DMRichPlot(RegionEnrichment, "genic", FALSE)
pdf("Gene_Region_Plot.pdf", width=12)
Region_plot
dev.off()

## CpG Annotation Enrichment Test ----------------------------------------------
CpG_annotation <- DMRichCpG(sigStats, statsbed, "hg19")

CpG_plot <- DMRichPlot(CpG_annotation, "CpG", FALSE)
pdf("CpG_Annotation_Plot.pdf", width=12)
CpG_plot
dev.off()

## Combined Enrichments Plots ----------------------------------------

# figure out way to add Chromatin marks to plot after chromHMM

CpG_annotation <- as.data.frame(CpG_annotation)
# write_xlsx(CpG_annotation, "./DMRichments/All DMRs_CpG_enrichments.xlsx")

RegionEnrichment <- as.data.frame(RegionEnrichment)
# write_xlsx(RegionEnrichment, "./DMRichments/All DMRs_genic_enrichments.xlsx")

# from .tsv output of chromHMM in DMRichR_ChromHMM.R
ChromHMMEnrichment <- read.table(file = 'allEnrichments_ChromHMM.tsv', sep = '\t', header = TRUE)

ChromHMMEnrichment <- ChromHMMEnrichment %>% dplyr::mutate(antibody = dplyr::recode_factor(antibody,
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
)

ChromHMMEnrichment <- ChromHMMEnrichment[,c('antibody', 'oddsRatio', 'qValue')]
colnames(ChromHMMEnrichment) <- c('Annotation', 'OR', 'fdr')

# combine both enrichments into single tibble object
CombinedEnrichment <- rbind(CpG_annotation, RegionEnrichment)
CombinedEnrichment <- CombinedEnrichment[, c('Annotation', 'OR', 'fdr')]
CombinedEnrichment <- rbind(CombinedEnrichment, ChromHMMEnrichment)
CombinedEnrichment$Direction <- "all"
CombinedEnrichment <- as_tibble(CombinedEnrichment)

# use the genic option to plot
CombinedPlot <- DMRichPlot(CombinedEnrichment, "genic", TRUE)
pdf("Combined_Annotation_Plot_chromHMM.pdf", width=12)
CombinedPlot
dev.off()


