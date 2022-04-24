# BiocManager::install("qusage")
library(limma)
library(qusage)
library(data.table)
library(plyr)
library(dplyr)
library(biomaRt)
library(missMethyl)
library(org.Hs.eg.db)
library(readxl)

# Update baseDir to proper location
baseDir <- "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN"

# set working directory to be folder containing root
setwd(baseDir)

################## For Methylation ####################
# load in differentially methylated probes
df_DMP <- read.csv("./results/data/DMPs_anno.csv")
df_DMP_sig <- filter(df_DMP, pval < 1e-7)

# load in most recent annotation file from illumina
# https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html
# v1.0 B4 manifest file same one used to annotate CpGs in missMethyl/minfi
minann <- read.csv("./data/datasets/MethylationEPIC_v-1-0_B4.csv", skip = 7)
rownames(minann) <- minann$IlmnID

# get annotation to join sig DMPs
#### NOTE #### run .getFlatAnnotation function in code_snippets.R before running below
ann <- .getFlatAnnotation("EPIC", anno = minann) # 1483 of 4242 CpGs have gene associated
# ann <- .getFlatAnnotation("EPIC")
ann$X <- rownames(ann)

# merge sig DMPs with mapping to Entrez Gene IDs
df <- merge(df_DMP_sig, ann, by="X")
df <- df %>% arrange(pval)
df$abs.tscore <- abs(df$t)
df$ranked.tvals <- order(df$abs.tscore, decreasing = T)

################# For RNA ######################
# load in differentially expressed genes
df_RNA <- read.csv("./results/data/EdgeR_RNA_sig_genes.csv")

# use biomaRt to find mapping info between ensembl and entrez gene IDS
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) # if unresponsive run: httr::set_config(httr::config(ssl_verifypeer = FALSE))
genes <- getBM(
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  mart=mart,
  useCache = FALSE)

df <- merge(df_RNA, genes, by.x="X", by.y="ensembl_gene_id")
df <- df %>% arrange(PValue)
df$abs.tscore <- abs(df$logFC)
# df$abs.tscore <- df$LR

# remove duplicate entrez gene id two ensembl match to same gene
df <- df[!df$X == "ENSG00000262102",]
df <- df[!is.na(df$entrezgene_id),] # 1 missing entrezgene_id
rownames(df) <- df$entrezgene_id
colnames(df)[colnames(df) == "entrezgene_id"] <- "entrezid"

df$ranked.tvals <- order(df$abs.tscore, decreasing = T)

################ For Both ###################

# load in datasets
# Transcription Factor
tft.sets <- read.gmt('./data/datasets/c3.tft.v7.5.1.entrez.gmt')

# Reactome
reactome.sets <- read.gmt('./data/datasets/c2.cp.reactome.v7.5.1.entrez.gmt')

# mitoCarta
mt.genes <- as.data.frame(read_excel("./data/datasets/Human.MitoCarta3.0.xlsx", sheet = "B Human All Genes"))
mt.pathways <- as.data.frame(read_excel("./data/datasets/Human.MitoCarta3.0.xlsx", sheet = "C MitoPathways"))

mt.sets <- as.list(strsplit(mt.pathways$Genes, split = ", "))
names(mt.sets) <- mt.pathways$MitoPathway
for (i in 1:length(mt.sets)) {
	set <- as.data.frame(unique(mt.sets[[i]]))
	colnames(set) <- c("set")
	set.entrez <- merge(set, mt.genes, by.x = "set", by.y = "Symbol")
	mt.sets[[i]] <- set.entrez$HumanGeneID
}

# perform.t.tests adapted from https://github.com/ArkingLab/mtDNA_GE_scripts/blob/master/Paper.Scripts/5.%20GO%20enrichment/permute.GO.all.tissues.R

perform.t.tests <- function(gene.sets, with.gene){
		signif.gene.sets <- as.data.frame(matrix(nrow = 1, ncol = 7))
		i <- 1
		colnames(signif.gene.sets) <- c('Gene.Set.Name', 'T.test.pval', 'Rank.t.test.pval', 'Num.genes.in.set', 'Beta', 'Confint.Upper', 'Confint.Lower')

		for(i in 1:length(gene.sets)){
			test.set <- gene.sets[[i]]
			set.name <- names(gene.sets[i])
			selected.indices <- which(with.gene$entrezid %in% test.set)
			in.set <- with.gene[selected.indices,]
			out.set <- with.gene[-selected.indices,]

			
			if(nrow(in.set) < 3){ # if this is set to 2, for the ranked.t.stats some of the results will be 0!!!!
				# print(paste0("Not enough samples for ", set.name))
			} else{
				t.stats <- t.test(in.set$abs.tscore, out.set$abs.tscore)
				ranked.t.stats <- t.test(in.set$ranked.tvals, out.set$ranked.tvals)
				sig <- c(set.name, t.stats$p.value, ranked.t.stats$p.value, nrow(in.set), t.stats$estimate[1]-t.stats$estimate[2], t.stats$conf.int[1], t.stats$conf.int[2])
				signif.gene.sets <- rbind(signif.gene.sets, sig)
				i <- i + 1
				if(i %% 100 == 0){
				 	print(paste0('On term ', i , ' out of ', length(gene.sets)))
				}
			}
		}

		signif.gene.sets <- na.omit(signif.gene.sets)
		signif.gene.sets$T.test.pval <- as.numeric(signif.gene.sets$T.test.pval)

		signif.gene.sets <- signif.gene.sets[order(signif.gene.sets$T.test.pval, decreasing = F),]
		return(signif.gene.sets)
	}

all.tft.sets <- perform.t.tests(tft.sets, df)
write.csv(all.tft.sets, './results/data/overrepresentation_tft_meth.csv')
all.reactome.sets <- perform.t.tests(reactome.sets, df)
write.csv(all.reactome.sets, './results/data/overrepresentation_reactome_meth.csv')
all.mt.sets <- perform.t.tests(mt.sets, df)
write.csv(all.mt.sets, './results/data/overrepresentation_mt_meth.csv')

with.gene.permute <- df
set.seed(1)

#### Get permutation cutoff for KEGG ####

start <- Sys.time()
all.min.tft <- numeric()
all.min.reactome <- numeric()
all.min.mt <- numeric()

for(p in 1:100){
  with.gene.permute$abs.tscore <- sample(with.gene.permute$abs.tscore)

  ### tft
  permute.tft <- perform.t.tests(tft.sets, with.gene.permute)
  min.pval <- permute.tft$T.test.pval[1]
  all.min.tft <- c(all.min.tft, min.pval)

  ### reactome
  permute.reactome <- perform.t.tests(reactome.sets, with.gene.permute)
  min.pval <- permute.reactome$T.test.pval[1]
  all.min.reactome <- c(all.min.reactome, min.pval)

  ### mt
  permute.mt <- perform.t.tests(mt.sets, with.gene.permute)
  min.pval <- permute.mt$T.test.pval[1]
  all.min.mt <- c(all.min.mt, min.pval)

  print(paste0('On permutation ', p))
}
end <- Sys.time()
runtime <- end-start

print(paste0('Total runtime: ', runtime))

### save permuted values
all.perm.cutoff.forset <- data.frame(TFT = all.min.tft, REACTOME = all.min.reactome, MT = all.min.mt)

all.perm.cutoff.forset$TFT <- all.perm.cutoff.forset$TFT[order(all.perm.cutoff.forset$TFT, decreasing = F)]
all.perm.cutoff.forset$REACTOME <- all.perm.cutoff.forset$REACTOME[order(all.perm.cutoff.forset$REACTOME, decreasing = F)]
all.perm.cutoff.forset$MT <- all.perm.cutoff.forset$MT[order(all.perm.cutoff.forset$MT, decreasing = F)]

# save(all.perm.cutoff.forset, file = 'all.perm.cutoff.forset.rds')

sig.tft.results <- subset(all.tft.sets, T.test.pval < all.perm.cutoff.forset$TFT[nrow(all.perm.cutoff.forset) * 0.05])	
sig.reactome.results <- subset(all.reactome.sets, T.test.pval < all.perm.cutoff.forset$REACTOME[nrow(all.perm.cutoff.forset) * 0.05])	
sig.mt.results <- subset(all.mt.sets, T.test.pval < all.perm.cutoff.forset$MT[nrow(all.perm.cutoff.forset) * 0.05])

if(nrow(sig.tft.results) != 0){sig.tft.results$Set <- 'TFT'}		
if(nrow(sig.reactome.results) != 0){sig.reactome.results$Set <- 'REACTOME'}	
if(nrow(sig.mt.results) != 0){sig.mt.results$Set <- 'MT'}	

all.results <- plyr::rbind.fill(sig.tft.results, sig.reactome.results)
# save(all.results, file = 'sig.pathways.nopseudo.permutecut.rds')
write.csv(all.results, "./results/data/tft_reactome_mito_meth.csv")

############################### mitoCarta genes ##############################
# Find list of significant genes that are in mitoCarta

# read in RNA sig genes
df_RNA <- read.csv("./results/data/EdgeR_RNA_sig_genes.csv")

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) # if unresponsive run: httr::set_config(httr::config(ssl_verifypeer = FALSE))
genes <- getBM(
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  mart=mart,
  useCache = FALSE)

df_RNA <- merge(df_RNA, genes, by.x="X", by.y="ensembl_gene_id")

mt.genes <- as.data.frame(read_excel("./data/datasets/Human.MitoCarta3.0.xlsx", sheet = "A Human MitoCarta3.0"))

# get data for output
sig_genes <- df_RNA$X
sig_genes <- sig_genes[!is.na(sig_genes)]
mt_genes <- mt.genes$EnsemblGeneID_mapping_version_20200130
mt_genes <- mt_genes[!is.na(mt_genes)]
mt_gene_symbol <- mt.genes$Symbol
mt_gene_symbol <- mt_gene_symbol[!is.na(mt_gene_symbol)]
mt_gene_desc <- mt.genes$Description
mt_gene_desc <- mt_gene_desc[!is.na(mt_gene_desc)]

entrez <- character(length(sig_genes))
symbol <- character(length(sig_genes))
desc <- character(length(sig_genes))
for (i in 1:length(sig_genes)){
    for (j in 1:length(mt_genes)){
        if (sig_genes[i] == mt_genes[j]) {
            entrez[i] = mt_genes[j]
			symbol[i] = mt_gene_symbol[j]
			desc[i] = mt_gene_desc[j]
        } 
    }
}

entrez <- entrez[!entrez == ""]
symbol <- symbol[!symbol == ""]
desc <- desc[!desc == ""]
df <- data.frame(entrez, symbol, desc)
colnames(df) <- c("Entrez_GeneID", "Symbol", "Description")
write.csv(df, "./results/data/mitoCarta_sig_genes.csv")

# 5 of 179 significant genes in mitocarta
# 1002 of 14149 genes in mitocarta
df_RNA_All <- read.csv("./results/data/EdgeR_RNA_all_genes.csv")
df_RNA_All <- merge(df_RNA_All, genes, by.x="X", by.y="ensembl_gene_id")

sig_genes_all <- df_RNA_All$entrezgene_id
sig_genes_all <- sig_genes_all[!is.na(sig_genes_all)]

entrez_all <- logical(length(sig_genes_all))
for (i in 1:length(sig_genes_all)){
    for (j in 1:length(mt_genes)){
        if (sig_genes_all[i] == mt_genes[j]) {
            entrez_all[i] = TRUE
        } 
    }
}

library(stats)
M <- as.table(rbind(c(5, 1002), c(174, 13147)))
dimnames(M) <- list(mitoCarta = c("T", "F"), count = c("Significant", "Not"))

chi <- chisq.test(M)