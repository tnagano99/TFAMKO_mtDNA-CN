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

############# Transcription Factor Fisher combined ###############################
TFT <- read.csv("TFT_Probe_RNA_Fisher.csv")
TFT <- TFT[, c("Gene.Set.Name", "P.DE_Meth", "P.DE_RNA", "Fisher")]
TFT$P.DE_Meth <- format(as.numeric(TFT$P.DE_Meth), digits=3, scientific=T)
TFT$P.DE_RNA <- format(as.numeric(TFT$P.DE_RNA), digits=3, scientific=T)
TFT$Fisher <- format(as.numeric(TFT$Fisher), digits=3, scientific=T)
kable(TFT[1:10,], caption="MetricsCompare", digits=4) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width=F) %>% save_kable("TopTFT_DMP_RNA_Fisher.png")

################## DMRichR and chromHMM results ##################################
Comb_Enrich <- read.csv("./DMRichR/Functional_Enrichments_DMRichR.csv")
Comb_Enrich <- Comb_Enrich %>% arrange(p.value)
Comb_Enrich$X <- NULL
Comb_Enrich$Direction <- NULL
Comb_Enrich$OR <- format(as.numeric(Comb_Enrich$OR), digits=3, scientific=T)
Comb_Enrich$CIlower <- format(as.numeric(Comb_Enrich$CIlower), digits=3, scientific=T)
Comb_Enrich$CIupper <- format(as.numeric(Comb_Enrich$CIupper), digits=3, scientific=T)
Comb_Enrich$p.value <- format(as.numeric(Comb_Enrich$p.value), digits=3, scientific=T)
Comb_Enrich$fdr <- format(as.numeric(Comb_Enrich$fdr), digits=3, scientific=T)
kable(Comb_Enrich, caption="MetricsCompare", digits=4) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width=F) %>% save_kable("TopFunctional_DMRichR.png")

#################### Transcription Factor, Reactome and MitoCarta ###########################
overRep <- read.csv("overrepresentation_mt_log.csv")
overRep <- overRep[,c("Gene.Set.Name", "T.test.pval", "Num.genes.in.set", "Beta", "Confint.Upper", "Confint.Lower")]
overRep <- overRep %>% arrange(T.test.pval)
overRep$T.test.pval <- format(as.numeric(overRep$T.test.pval), digits=3, scientific=T)
overRep$Beta <- format(as.numeric(overRep$Beta), digits=3, scientific=F)
overRep$Confint.Upper <- format(as.numeric(overRep$Confint.Upper), digits=3, scientific=F)
overRep$Confint.Lower <- format(as.numeric(overRep$Confint.Lower), digits=3, scientific=F)
kable(overRep[1:10,], caption="MetricsCompare", digits=4) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width=F) %>% save_kable("TopOverRep_mt_RNA.png")