# TFAMKO_mtDNA-CN

Repository for MHI 4980E Thesis Project as part of the Castellani Lab

**Topic:** Effect of CRISPR Induced Mitochondrial DNA Varition on the Nuclear Epigenome and Transcriptome

## Table of Contents

- Overview
- Hypothesis
- Objectives
- Folder Structure
- File Descriptions

### Overview

The main objective of this project is to identify underlying biological mechanisms where variation in mitochondrial DNA copy number (mtDNA-CN) affect nuclear DNA (nDNA) methylation in a mitochondrial transcription factor A (TFAM) heterozygous knockout cell line model. The epigenome and transcriptome data generated from the cell line model will be analyzed and integrated, using bioinformatics programs on the command line and R, to determine sites of differential methylation and differential gene expression to determine relevant biological mechanisms.

#### Hypothesis

In vitro reduction of mtDNA-CN leads to differentially methylated sites/regions and differential gene expression in shared genes and biological pathways which mediate the effect of mtDNA-CN on disease

#### Objectives

- Analyze data generated from cell model to determine differentially methylated sites and differential gene expression. 
- Integrate results of methylation and gene expression analyses to find overlapping pathways and genes
- Perform functional enrichment analyses on significant results to discover relevant biological mechanisms
    - Functional Enrichment Analyses
    - Feature Overrepresentation

### Folder Structure

This project has four main folders:
- data: contains all data generated from cell line, annotation files, and database files
- doc: contains documents related to thesis
- results: contains all outputs of scripts in tables (data) or visualizations (plots)
- src: contains all R scripts and shell scripts for all analyses and visualizations

**Subfolders**

data
- datasets
- fasta
- Roadmap_ChromHMM
- roadmap_epigenomics
- run_combined

doc

results
- data
    - CpG_mtDNA-CN
    - DMRcate
    - DMRichR
    - EdgeR
    - ELMER
    - manhattan
    - minfi
    - pathview
    - QQ_plots
- plots
    - CpG_mtDNA-CN
    - DMRcate
    - DMRichR
    - EdgeR
    - ELMER
    - manhattan
    - minfi
    - pathview
    - QQ_plots
src

### File Descriptions

Below are short descriptions of the function of each script in the src folder

**EPIC_methylation_minfi_combined_runs**

Note: EPIC_methylation_minfi.R was previously used when running normalization and preprocessing for the two runs separately. no longer used
1. minfi 
    - reads in intensity values from .idat files taken from the Illumina Infinium EPIC BeadChip 
    - performs normalization and quality control by removing poor performing probes
    - plots to visualize normalization at end of script
    - plots saved in plots/minfi
    - creates **beta.csv** and **mVal.csv**
2. DMPFinder
    - determines differentially methylated sites
    - generates summary stats saved in **dmp_cont.csv**
3. DMRcate
    - determines differentially methylated regions
    - generates summary stats saved in **DMRS_anno.csv**
    - functions to visualize regions at end of script
    - plots saved in plots/DMRcate

**EdgeR**

1. Convert Kallisto Output
    - take outputs from Kallisto and convert to gene-level counts
2. EdgeR
    - performs quality control and normalization for library size
    - determines differentially expressed genes
    - generates summary stats saved in **EdgeR_RNA_sig_genes.csv** and **EdgeR_RNA_all_genes.csv**
3. tximport
    - standard package to convert Kallisto output to gene-level counts

**ELMER-Analysis**

Note the ELMER-Analysis-DMP-Cont.R was used in ELMER.sh to run the script as a job on Compute Canada. 
1. ELMER
    - takes in DNA methylation and gene expression data
        - use significant results from previous analyses
    - finds nearest 20 genes to each differentially methylated sites
    - tests for inverse correlations between DNA methylation and gene expression
    - significant gene-probe pairs saved in **getPair.ALL_DMP_CONT_EDGER_MAX_SIG.pairs.significant.filtered.csv**
2. Visualizations
    - many different options to visualize relationships
    - all results in plots/ELMER

**Meth_Analysis_missMethyl**

1. missMethyl
    - perform functional enrichments for differentially methylated sites and differentially methylated regions
    - tests for enrichments in the GO (Gene Ontology) and KEGG (Kyoto Encyclopedia of Genes and Genomes)
    - differentially methylated sites results saved as **GO_Cont.csv** and **KEGG_Cont.csv**
    - differentially methylated regions results saved as **GO_DMRcate.csv** and **KEGG_DMRcate.csv**

**RNA_Analysis**

1. limma
    - perform functional enrichments for differentially expressed genes
    - tests for enrichments in the GO (Gene Ontology) and KEGG (Kyoto Encyclopedia of Genes and Genomes)
    - results saved as **GO_RNA_LRT_EDGER.csv** and **KEGG_RNA_LRT_EDGER.csv**

**Kegg_pathview**

1. pathview
    - creates annotated visualizations of KEGG pathways with amount of methylation or gene expression between knockout and control
    - creates all pathview plots; plots named according to KEGG identifier

**Feature_Overepresentation**

1. Transcription Factor, Reactome, MitoCarta
    - perform functional enrichment analyses for differentially methylated sites, differentially methylated regions and differentially expressed genes in the Transcription Factor and Reactome databases taken from MSigDB and the mitoCarta 3.0 database
    - outputs all results with the name overrepresentation*.csv

**Fisher_Combined.R**

1. metap
    - use the metap function to perform Fisher's Combined test to meta-analyze the methylation and RNA functional enrichment results together
    - performed for GO, KEGG, Transcription Factor, Reactome and mitoCarta

**DMRich and DMRich_ChromHMM**

1. DMRich
    - perform functional enrichments of functional groups using the methylation data
    - checks for enrichment for CpGs and for functional gene groups

2. ChromHMM and Roadmap Epigenomics
    - perform functional enrichments of epigenomic elements using the methylation data
    - checks for enrichment in 15-state chromatin state model and the related 5 core histone modifications

**format_tables**

1. Format all tables
    - creates nicely formatted tables for all outputs included in powerpoints in doc/thesis_presentation.ppt

**code_snippets**

1. manhattanraw
    - creates manhattan plots from methylation and RNA data
    - outputs saved in plots/manhattan
2. QQ Plot and Lambda Calculation
    - creates qq-plots and calculates lambda value from p-values of summary stats for all analyses
    - outputs saved in plots/QQ_plots
3. cpg.annotate
    - modified function from DMRcate to use p-value as cutoff instead of fdr
    - run this block before cpg.annotate for DMRcate analysis
4. getFlatAnnotation
    - taken from missMethyl source code to read in annotation file in Feature_Overepresentation.R
5. Linear Mixed Model Code:
    - **Not used**
    - runs a linear mixed model on the methylation data

**plots**

1. Manhattan Plots
    - uses the manhattanraw function from code_snippets.R to generate manhattan plots

2. QQ Plots
    - uses the qq function from qqman to generate qq-plots

**No Longer Used**

Below are the list of files not used in the analysis. You can ignore all of these files

- get_meme_sequences.R
- meme_commands.txt
- Meme.sh
- run_linear_model.sh
- ELMER-Analysis-DMP-Cont.R
- EPIC_methylation_minfi.R
- linear_model_mtDNACN.R
     