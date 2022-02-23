library(magrittr)
library(kableExtra)
library("magick")
library("webshot")

setwd("Thesis")

GO <- read.csv("GO_DMRcate_pval.csv")

GO$P.DE <- format(as.numeric(GO$P.DE), digits=3, scientific=T)
# GO$GO_label <- NULL
GO$FDR <- format(as.numeric(GO$FDR), digits=3, scientific=T)
kable(GO[1:10,], caption="MetricsCompare", digits=4) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width=F) %>% save_kable("TopGO_DMRcate_new.png")

KEGG <- read.csv("KEGG_DMRcate_pval.csv")
KEGG$P.DE <- format(as.numeric(KEGG$P.DE), digits=3, scientific=T)
# KEGG$KEGG_label <- NULL
KEGG$FDR <- format(as.numeric(KEGG$FDR), digits=3, scientific=T)
kable(KEGG[1:6,], caption="MetricsCompare", digits=4) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width=F) %>% save_kable("TopKEGG_DMRcate_new.png")

DMPs <- read.csv("DMPs_anno.csv")
DMPs <- DMPs[, c("X", "intercept", "beta", "t", "pval", "qval", "chr", "pos", "GencodeCompV12_NAME", "GencodeCompV12_Accession", "GencodeCompV12_Group")]

DMPs$pval <- format(as.numeric(DMPs$pval), digits=3, scientific=T)
DMPs$qval <- format(as.numeric(DMPs$qval), digits=3, scientific=T)
kable(DMPs[1:15,], caption="MetricsCompare", digits=4) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width=F) %>% save_kable("TopDMPs.png")

DMRs <- read.csv("DMRS_anno.csv")

DMRs$Fisher <- format(as.numeric(DMRs$Fisher), digits=3, scientific=T)
kable(DMRs[1:15,], caption="MetricsCompare", digits=4) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width=F) %>% save_kable("TopDMRs.png")

RNAseq <- read.csv("SleuthAllGenesAnnotatedRNASeqResultsGeneWise_cleaned.csv")

RNAseq <- RNAseq[, c("target_id","test_stat", "pval", "ext_gene", "chr", "description")]
RNAseq$pval <- format(as.numeric(RNAseq$pval), digits=3, scientific=T)
kable(RNAseq[1:15,], caption="MetricsCompare", digits=4) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width=F) %>% save_kable("TopRNA.png")
