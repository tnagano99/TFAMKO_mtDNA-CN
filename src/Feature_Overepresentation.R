# BiocManager::install("qusage")
library(limma)
library(qusage)
library(data.table)
library(plyr)
library(dplyr)
library(biomaRt)
library(missMethyl)
library(org.Hs.eg.db)

# Update baseDir to proper location
baseDir <- "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN"

# set working directory to be folder containing root
setwd(baseDir)

# load in differentially expressed genes
# df_RNA <- read.csv("./results/data/EdgeR_RNA_all_genes.csv")

# load in differentially methylated probes
df_DMP <- read.csv("./results/data/DMPs_anno.csv")
df_DMP_sig <- filter(df_DMP, pval < 1e-7)

# load in most recent annotation file from illumina
# https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html
# v1.0 B5 manifest file
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

# get CpGs mapped to Entrez Gene IDs using missMethyl package
# df_DMP_sig_genes <- getMappedEntrezIDs(df_DMP_sig$X, array.type = "EPIC")

# load in datasets
tft.sets <- read.gmt('./data/datasets/c3.tft.v7.5.1.entrez.gmt')
reactome.sets <- read.gmt('./data/datasets/c2.cp.reactome.v7.5.1.entrez.gmt')

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
# write.csv(all.tft.sets, './results/data/overrepresentation_tft.csv')
all.reactome.sets <- perform.t.tests(reactome.sets, df)
# write.csv(all.reactome.sets, './results/data/overrepresentation_reactome.csv')

with.gene.permute <- df
set.seed(1)

#### Get permutation cutoff for KEGG ####

start <- Sys.time()
all.min.tft <- numeric()
all.min.reactome <- numeric()

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

  print(paste0('On permutation ', p))
}
end <- Sys.time()
runtime <- end-start

print(paste0('Total runtime: ', runtime))

### save permuted values
all.perm.cutoff.forset <- data.frame(TFT = all.min.tft, REACTOME = all.min.reactome)

all.perm.cutoff.forset$TFT <- all.perm.cutoff.forset$TFT[order(all.perm.cutoff.forset$TFT, decreasing = F)]
all.perm.cutoff.forset$REACTOME <- all.perm.cutoff.forset$REACTOME[order(all.perm.cutoff.forset$REACTOME, decreasing = F)]

# save(all.perm.cutoff.forset, file = 'all.perm.cutoff.forset.rds')

sig.tft.results <- subset(all.tft.sets, T.test.pval < all.perm.cutoff.forset$TFT[nrow(all.perm.cutoff.forset) * 0.05])	
sig.reactome.results <- subset(all.reactome.sets, T.test.pval < all.perm.cutoff.forset$REACTOME[nrow(all.perm.cutoff.forset) * 0.05])	


if(nrow(sig.tft.results) != 0){sig.tft.results$Set <- 'TFT'}		
if(nrow(sig.reactome.results) != 0){sig.reactome.results$Set <- 'REACTOME'}	

all.results <- plyr::rbind.fill(sig.tft.results, sig.reactome.results)
# save(all.results, file = 'sig.pathways.nopseudo.permutecut.rds')