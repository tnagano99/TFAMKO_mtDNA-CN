library(magrittr)
library(kableExtra)
library("magick")
library("webshot")

# set working directory to folder with all results
setwd("results/data")

################ standard table with p-values and fdr #######################
# read in the GO Result update below for DMP, DMR, RNA
GO <- read.csv("GO_RNA_LRT_EDGER.csv")

GO$P.DE <- format(as.numeric(GO$P.DE), digits=3, scientific=T)
GO$GO_label <- NULL # run for RNA
# GO$FDR <- format(as.numeric(GO$FDR), digits=3, scientific=T) # run for DMP and DMR
kable(GO[1:10,], caption="MetricsCompare", digits=4) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width=F) %>% save_kable("TopGO_RNA_EDGER.png")

# read in the KEGG Results; update below for DMP, DMR, RNA
KEGG <- read.csv("KEGG_RNA_LRT_EDGER.csv")
KEGG$P.DE <- format(as.numeric(KEGG$P.DE), digits=3, scientific=T)
KEGG$KEGG_label <- NULL # run for RNA
# KEGG$FDR <- format(as.numeric(KEGG$FDR), digits=3, scientific=T) # run for DMP and DMR
kable(KEGG[1:6,], caption="MetricsCompare", digits=4) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width=F) %>% save_kable("TopKEGG_RNA_EDGER.png")

################ raw DMPs and DMRs with annotation ######################

DMPs <- read.csv("DMPs_anno.csv")
DMPs <- DMPs[, c("X", "intercept", "beta", "t", "pval", "qval", "chr", "pos", "GencodeCompV12_NAME", "GencodeCompV12_Accession", "GencodeCompV12_Group")]
DMPs$pval <- format(as.numeric(DMPs$pval), digits=3, scientific=T)
DMPs$qval <- format(as.numeric(DMPs$qval), digits=3, scientific=T)
kable(DMPs[1:15,], caption="MetricsCompare", digits=4) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width=F) %>% save_kable("TopDMPs.png")

DMRs <- read.csv("DMRS_anno.csv")
DMRs$Fisher <- format(as.numeric(DMRs$Fisher), digits=3, scientific=T)
kable(DMRs[1:15,], caption="MetricsCompare", digits=4) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width=F) %>% save_kable("TopDMRs.png")

################### RNA significant genes ###########################

RNAseq <- read.csv("EdgeR_RNA_sig_genes.csv")
# RNAseq <- RNAseq[, c("target_id","test_stat", "pval", "ext_gene", "chr", "description")] # for sleuth results
# RNAseq <- RNAseq[, c("X","LR", "PValue", "FDR")] # for sleuth results
RNAseq$PValue <- format(as.numeric(RNAseq$PValue), digits=3, scientific=T)
RNAseq$FDR <- format(as.numeric(RNAseq$FDR), digits=3, scientific=T)
kable(RNAseq[1:15,], caption="MetricsCompare", digits=4) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width=F) %>% save_kable("TopRNAEdgeR.png")

########## GO Fisher combined results ###################
GO <- read.csv("GO_Probe_RNA_Fisher_EDGER.csv") # change between Probe and Region
GO <- GO[, c("GO_label", "TERM", "P.DE_Meth", "P.DE_RNA", "Fisher")]

GO$P.DE_Meth <- format(as.numeric(GO$P.DE_Meth), digits=3, scientific=T)
GO$P.DE_RNA <- format(as.numeric(GO$P.DE_RNA), digits=3, scientific=T)
GO$Fisher <- format(as.numeric(GO$Fisher), digits=3, scientific=T)
kable(GO[1:10,], caption="MetricsCompare", digits=4) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width=F) %>% save_kable("TopGO_DMR_RNA_Fisher_EDGER.png")

########## KEGG Fisher combined results #################
KEGG <- read.csv("KEGG_Region_RNA_Fisher_EDGER.csv") # change between Probe and Region
KEGG <- KEGG[, c("KEGG_label", "Description", "P.DE_Meth", "P.DE_RNA", "Fisher")]
KEGG$P.DE_Meth <- format(as.numeric(KEGG$P.DE_Meth), digits=3, scientific=T)
KEGG$P.DE_RNA <- format(as.numeric(KEGG$P.DE_RNA), digits=3, scientific=T)
KEGG$Fisher <- format(as.numeric(KEGG$Fisher), digits=3, scientific=T)
kable(KEGG[1:10,], caption="MetricsCompare", digits=4) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width=F) %>% save_kable("TopKEGG_DMR_RNA_Fisher_EDGER.png")
