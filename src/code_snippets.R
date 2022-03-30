# Code snippets in this file belong to Dr. Castellani of the Castellani Lab at Schulich School of Medicine and Dentistry at Western University

# Code to Create Manhattan plots

EPIC_CC <- read.table("/dcs01/arking/arkinglab/active/projects/aric/epigenetics/ccastell/Whites/EPICResults2.txt", header=T)
setwd("/dcs01/arking/arkinglab/active/projects/aric/epigenetics/ccastell/EPIC/")
annot <- read.csv("MethylationEPIC_15073387_v-1-0.csv", skip=7)

mergeManhattan <- merge(EPIC_CC, annot, by.x="row.names", by.y="IlmnID")
DMS <- mergeManhattan

colnames(DMS)[2] <- "estimate"
colnames(DMS)[5] <- "P.Value"
colnames(DMS)[16] <- "chr"
colnames(DMS)[17] <- "start"


manhattanraw<-function(DMS, filename, sig=NULL){
	# DMS needs to be a dataframe with columns 
  # P.Value: the p-values of the differentially methylated sites
  # chr: the chromosome the CpG site is on
  # start: the position on the chromosome of the CpG site

	#if necessary, get rid of sex chromosomes
	data<-DMS
	data<-subset(data, (data$chr!="0"))
	data<-subset(data, (data$start!="NA"))
	data<-subset(data, (data$chr!="X"))
	data<-subset(data, (data$chr!="Y"))

	data$chr <- gsub("chr","",data$chr)
	data$chr <- as.numeric(data$chr)
	data2 <- data[order(data$chr),]
	data=data2

	# set the max value needed for the Manhattan plot
	ymin1= round( max(-log10(data$P.Value))+1)

  # set the title for the plot
	title=c()

  # initialize the jpeg file to create the plot
	jpeg(paste("./results/plots/manhattan/manhattan_", filename, ".jpeg",sep=""),res=400,width = 40, height = 12,units="cm")
	
  # initialize x-axis scale
  chr <- c(1:22)

	# set the position using start
	data$position<-round(data$start,digits=0)
	print(table(data$chr))

  # set the margins for the plot
	par(mar=c(5,5,2,2))

  # find the cutoff points for chromosomes
	phy.max<-tapply(data$start, data$chr,max,na.rm=T)
	cumlen=0

  # remove values without chromsome specified
	data <- data[!is.na(data$chr),]

  # determine actual position for CpGs using median
	for(i in chr){
		cat(paste("Now working on chromosome ",i,"\r"))
		data[data$chr==i,"loc"]<-data[data$chr==i,"position"]+cumlen
		cumlen<-cumlen+phy.max[i]
	}
	phy.med<-tapply(data$loc,data$chr,median,na.rm=T)[chr]
	print(phy.med)

  # convert p-values to -log10 transform
	data$mlgpval<- -log(data[,"P.Value"], base=10)

  # specify the plot parameters loc is the chromosome and mlgpval is the log10 p-value
	plot(data[,"loc"], data[,"mlgpval"],
    type="n",yaxt="n",xaxt="n",
    main=title,
		xlab="chromosome",
		ylab=expression(-log[10]*P),
		xlim=c(0,max(data$loc + 1,na.rm=T)),cex.lab=1.5,ylim=c(0,ymin1)
		col=ifelse(data$mlgpval< sig, (rep(c("black","gray48"),13)), "red" ))

  # set y axis labels and scale
	axis(side=2, at=seq(from=0,to=ymin1,by=2), labels=seq(from=0,to=ymin1,by=2), tick=T, cex.axis=0.9, las=1)

  # set x axis labels and scale
	axis(side=1, at=phy.med[c(1:22)], labels=c(1:22), cex.axis=0.8)

  # fill margins with white rectangles
	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
		col = "white")

  # output to console for updates
	for(i in chr){
		cat(paste("Now working on chromosome ",i,"\r"))
		points(data[data$chr==i,"loc"],data[data$chr==i,"mlgpval"],
			col=ifelse(data[data$chr==i,"mlgpval"] > sig, "red", col[i]),pch=20)
	}

	# add in line for significance
	if (is.null(sig)==FALSE){
		abline(h=sig,lty="dotted",lwd=2,col="chartreuse4")
	}

	dev.off()
}

######################### QQ Plot code ###################################

qq.chisq(summaryStats$pval)
# pvals <- summaryStats$pval
pvals <- summaryStats$PValue
observed <- sort(pvals)
observed2 <- c(length(pvals))
observed2null <- -(log10(observed2 / (length(observed2)+1)))
pvals <- c(pvals, observed2null)
observed <- sort(pvals)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
# creating uniform distribution
m="qqplot"
pdf(paste0("./results/plots/QQ_plots/QQ_EDGER.pdf"))
par(mfrow=c(1,1))
plot(c(0,20), c(0,20), col="red", lwd=4, type="l", xlab="Expected (-logP)",
	ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l", main=m)
points(lexp, lobs, pch=23, cex=.5, col="black", bg="black")
dev.off()

################ code to check lambda value of p-values ###################
library(QCEWAS)
P_lambda(pvalsModelA) # replace with vector of p-values

######################### DMRcate cpg.annotate ###########################################
# updated to use the p-value instead of FDR value
# run the block of code from lines 137 to 277 for DMR analysis before calling cpg.annotate in EPIC_methylation_minfi_combined_runs.R

cpg.annotate <- function (datatype = c("array", "sequencing"), object, what = c("Beta", "M"), 
          arraytype = c("EPIC", "450K"), analysis.type = c("differential", "variability", "ANOVA", "diffVar"),
          design, contrasts = FALSE, cont.matrix = NULL, fdr = 0.05, coef, varFitcoef=NULL, topVarcoef=NULL, ...) 
{
  analysis.type <- match.arg(analysis.type)
  what <- match.arg(what)
  arraytype <- match.arg(arraytype)
  if (datatype == "array") {
    stopifnot(class(object)[1] %in% c("matrix", "GenomicRatioSet"))
    if (is(object, "matrix")) {
      if (arraytype == "450K") {
        grset <- makeGenomicRatioSetFromMatrix(mat = object, 
                                               array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19", 
                                               what = what)
      }
      if (arraytype == "EPIC") {
        grset <- makeGenomicRatioSetFromMatrix(mat = object, 
                                               array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19", 
                                               what = what)
      }
    }
    else {
      grset <- object
    }
    object <- getM(grset)
    switch(analysis.type, differential = {
      stopifnot(is.matrix(design))
      if (!contrasts) {
        stopifnot(colnames(design)[1] == "(Intercept)")
      } else {
        stopifnot(!is.null(cont.matrix))
      }
      fit <- lmFit(object, design, ...)
      if (contrasts) {
        stopifnot(coef %in% colnames(cont.matrix))
        fit <- contrasts.fit(fit, cont.matrix)
      }
      fit <- eBayes(fit)
      tt <- topTable(fit, coef = coef, number = nrow(object))
      nsig <- sum(tt$P.Value < fdr)
      if (nsig == 0) {
        message("Your contrast returned no individually significant probes. Try increasing the fdr. Alternatively, set pcutoff manually in dmrcate() to return DMRs, but be warned there is an increased risk of Type I errors.")
      }
      if (nsig > 0 & nsig <= 100) {
        message(paste("Your contrast returned", nsig, 
                      "individually significant probes; a small but real effect. Consider manually setting the value of pcutoff to return more DMRs, but be warned that doing this increases the risk of Type I errors."))
      }
      if (nsig > 100) {
        message(paste("Your contrast returned", nsig, 
                      "individually significant probes. We recommend the default setting of pcutoff in dmrcate()."))
      }
      betafit <- lmFit(ilogit2(object), design, ...)
      if (contrasts) {
        betafit <- contrasts.fit(betafit, cont.matrix)
      }
      betafit <- eBayes(betafit)
      betatt <- topTable(betafit, coef = coef, number = nrow(object))
      m <- match(rownames(tt), rownames(betatt))
      tt$diff <- betatt$logFC[m]
      m <- match(rownames(object), rownames(tt))
      tt <- tt[m, ]
      anno <- getAnnotation(grset)
      stat <- tt$t
      annotated <- GRanges(as.character(anno$chr), IRanges(anno$pos, anno$pos), stat = stat,
                           diff = tt$diff, ind.fdr = tt$adj.P.Val, is.sig = tt$P.Value < fdr)
      names(annotated) <- rownames(tt)
    }, variability = {
      RSanno <- getAnnotation(grset)
      wholevar <- var(object)
      weights <- apply(object, 1, var)
      weights <- weights/mean(weights)
      annotated <- GRanges(as.character(RSanno$chr), IRanges(RSanno$pos, RSanno$pos), stat = weights,
                           diff = rep(0, nrow(object)), ind.fdr = rep(0, nrow(object)), is.sig = weights > quantile(weights, 0.95))
      names(annotated) <- rownames(object)
    }, ANOVA = {
      message("You are annotating in ANOVA mode: consider making the value of fdr quite small, e.g. 0.001")
      stopifnot(is.matrix(design))
      fit <- lmFit(object, design, ...)
      fit <- eBayes(fit)
      sqrtFs <- sqrt(fit$F)
      sqrtfdrs <- p.adjust(fit$F.p.value, method = "BH")
      nsig <- sum(sqrtfdrs < fdr)
      if (nsig == 0) {
        message("Your design returned no individually significant probes for ANOVA. Try increasing the fdr. Alternatively, set pcutoff manually in dmrcate() to return DMRs, but be warned there is an increased risk of Type I errors.")
      }
      if (nsig > 0 & nsig <= 100) {
        message(paste("Your design returned", nsig, 
                      "individually significant probes for ANOVA; a small but real effect. Consider manually setting the value of pcutoff to return more DMRs, but be warned that doing this increases the risk of Type I errors."))
      }
      if (nsig > 100) {
        message(paste("Your design returned", nsig, 
                      "individually significant probes for ANOVA. We recommend the default setting of pcutoff in dmrcate(). Large numbers (e.g. > 100000) may warrant a smaller value of the argument passed to fdr"))
      }
      anno <- getAnnotation(grset)
      stat <- sqrtFs
      annotated <- GRanges(as.character(anno$chr), IRanges(anno$pos, anno$pos), stat = stat,
                           diff = 0, ind.fdr = sqrtfdrs, is.sig = sqrtfdrs < fdr)
      names(annotated) <-  rownames(object)
    }, diffVar = {
      stopifnot(is.matrix(design))
      if (!contrasts) {
        stopifnot(colnames(design)[1] == "(Intercept)")
      } else {
        stopifnot(!is.null(cont.matrix))
      }
      fitvar <- varFit(object, design = design, coef=varFitcoef)
      if (contrasts) {
        fitvar <- contrasts.varFit(fitvar, cont.matrix)
      }
      tt <- topVar(fitvar, coef = topVarcoef, number = nrow(object))
      nsig <- sum(tt$Adj.P.Value < fdr)
      if (nsig == 0) {
        message("Your contrast returned no individually significant probes. Try increasing the fdr. Alternatively, set pcutoff manually in dmrcate() to return DVMRs, but be warned there is an increased risk of Type I errors.")
      }
      if (nsig > 0 & nsig <= 100) {
        message(paste("Your contrast returned", nsig, 
                      "individually significant probes; a small but real effect. Consider manually setting the value of pcutoff to return more DVMRs, but be warned that doing this increases the risk of Type I errors."))
      }
      if (nsig > 100) {
        message(paste("Your contrast returned", nsig, 
                      "individually significant probes. We recommend the default setting of pcutoff in dmrcate()."))
      }
      m <- match(rownames(object), rownames(tt))
      tt <- tt[m, ]
      anno <- getAnnotation(grset)
      stat <- tt$t
      annotated <- GRanges(as.character(anno$chr), IRanges(anno$pos, anno$pos), stat = stat,
                           diff = 0, ind.fdr = tt$Adj.P.Value, is.sig = tt$Adj.P.Value < fdr)
      names(annotated) <-  rownames(tt)  
    })
    annotated <- sort(annotated)
    return(new("CpGannotated", ranges=annotated))
    
  }
  if (datatype == "sequencing") {
    stop("Sequencing mode is deprecated for cpg.annotate(). Please use sequencing.annotate().")
  }
  else {
    message("Error: datatype must be one of 'array' or 'sequencing'")
  }
}

########################## getFlatAnnotation from MissMethyl ###########################

# Taken from https://github.com/Oshlack/missMethyl/blob/master/R/gometh.R

.getFlatAnnotation <- function(array.type=c("450K","EPIC"),anno=NULL)
  # flatten 450k or EPIC array annotation
  # Jovana Maksimovic
  # 18 September 2018
  # Updated 18 September 2018
  # Modified version of Belida Phipson's .flattenAnn code
{
  if(is.null(anno)){
    if(array.type=="450K"){
      anno <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
    } else {
      anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    }
  }
  
  # get rid of the non-CpG sites
  ann.keep<-anno[grepl("^cg",anno$Name),]
  
  # get rid of CpGs that are not annotated
  missing<-ann.keep$UCSC_RefGene_Name==""
  ann.keep<-ann.keep[!missing,]
  
  # get individual gene names for each CpG
  geneslist<-strsplit(ann.keep$UCSC_RefGene_Name,split=";")
  names(geneslist)<-rownames(ann.keep)
  
  grouplist<-strsplit(ann.keep$UCSC_RefGene_Group,split=";")
  names(grouplist)<-rownames(ann.keep)
  
  flat<-data.frame(symbol=unlist(geneslist),group=unlist(grouplist))
  flat$symbol<-as.character(flat$symbol)
  flat$group <- as.character(flat$group)
  
  flat$cpg<- substr(rownames(flat),1,10)
  
  #flat$cpg <- rownames(flat)
  flat$alias <- suppressWarnings(limma::alias2SymbolTable(flat$symbol))
  
  #eg <- toTable(org.Hs.egSYMBOL2EG)
  eg <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, 
                                               keys=AnnotationDbi::keys(org.Hs.eg.db), 
                                               columns=c("ENTREZID","SYMBOL"), 
                                               keytype="ENTREZID"))
  colnames(eg) <- c("gene_id","symbol")
  
  m <- match(flat$alias,eg$symbol)
  flat$entrezid <- eg$gene_id[m]
  flat <- flat[!is.na(flat$entrezid),]
  
  # keep unique cpg by gene name annotation
  id<-paste(flat$cpg,flat$entrezid,sep=".")
  d <- duplicated(id)
  flat.u <- flat[!d,]
  flat.u
  # This randomly samples only 1 gene ID for multimapping CpGs
  #.reduceMultiMap(flat.u)
}


########################## Linear Mixed Model code #####################################

########## NOTE: below code is not used in project #####################
# code to run a linear model for each probe

library(lme4)
## Define the lm function
runlme <- function(thisdat) {
   lme1 <- eval(parse(text=expression));    
   ##Get the summary of the model
   smodel = summary(lme1);
   return(smodel)
 }

varnames <- c("mtDNACN")
Bmat <- SEmat <- Tmat <- matrix(NA,nrow=nrow(EPICTFAMKO),ncol=1)
rownames(Bmat) <- rownames(SEmat) <- rownames(Tmat) <- rownames(EPICTFAMKO)
colnames(Bmat) <- paste("Estimate",varnames,sep=".")
colnames(SEmat) <- paste("Std.Error",varnames,sep=".")
colnames(Tmat) <- paste("t-value",varnames,sep=".")

for (i in 1:nrow(EPICTFAMKO)) { 
if (i %% 1000 == 0) {cat(paste("On probe ",i,"\n",sep=""))} #outputs every 1000 probes just to check in
thisExpr <- as.numeric(EPICTFAMKO[i,])
expression <- "lmer(thisExpr~mtDNACN + (1|Batch)+(1|ID), na.action=na.exclude, control = lmerControl(calc.derivs = FALSE), REML=FALSE)"

designmatrix <- data.frame(thisExpr, mtDNACN, Batch, ID)
lme1.out <- try(runlme(designmatrix),silent=F);

if (substr(lme1.out[1],1,5)!="Error") {
tabOut <- lme1.out$coefficients
Bmat[i,] <- tabOut[2,"Estimate"]
SEmat[i,] <- tabOut[2,"Std. Error"]
Tmat[i,] <- tabOut[2,"t value"]
  } else {
    cat('Error in LME of Probe',rownames(EPICTFAMKO)[i],"id",'\n')
    cat('Setting P-value=NA,Beta value=NA, and =NA\n');
    Bmat[i,] <- SEmat[i,] <- Tmat[i,] <- NA;
  }
}

warnings()

FinalResults <- cbind(Bmat, SEmat, Tmat)
FinalResults <- as.data.frame(FinalResults)
zscores <- FinalResults[,3]
pvalue <- pchisq(zscores**2,1,lower.tail=F)
min(pvalue)
FinalResults2 <- cbind(FinalResults,pvalue)

write.csv(FinalResults2,"EPIClmerResults.csv", quote=F)

###################################################################