library(biomaRt)
library(regioneR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)

# get all ensembl genes into granges object
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) # if unresponsive run: httr::set_config(httr::config(ssl_verifypeer = FALSE))

genes <- getBM(attributes=c("chromosome_name", "start_position", "end_position", "ensembl_gene_id"),
        mart = mart,
        values = c("protein_coding"))

ranges <- toGRanges(genes)

ranges <- keepSeqlevels(ranges, c(1:22,"X", "Y"), pruning.mode = "coarse")

# read in RNA-Seq analysis results to get significant genes
df <- read.csv("./results/data/SleuthAllGenesAnnotatedRNASeqResultsGeneWise_cleaned.csv", header=T) # Likelihood test results
names(df)[names(df) == "target_id"] <- "ensembl_gene_id"

# merge data with entrez gene ids
df_anno <- merge(df, genes, by = "ensembl_gene_id")
df_sig <- filter(df_anno, pval < 3.59e-6)

ranges <- ranges[(elementMetadata(ranges)[,1] %in% df_sig$ensembl_gene_id)]

flanking <- flank(ranges, 1000, start = TRUE)

seqlevelsStyle(flanking) <- "UCSC"

seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, flanking)

writeXStringSet(seqs, "/home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN/data/fasta/sig_genes_upstream_1000.fa", format="fasta")