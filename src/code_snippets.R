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


manhattan<-function(DMS, filename, sig=NULL){
	#require(methylKit)

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
	#data2$meth.diff <-
	#cex.val <- 0.8 - ((abs(data[,"mean.meth.diff"]) < 5 )*.5) +
	#((abs(data[,"mean.meth.diff"]) > 20)*.7)
	#data$cex.val<-cex.val

	# library(GenABEL)
	# z.sq<-(data$Effect / data$StdErr)^2
	# lambda<-estlambda(data$P.Value)
	# cat(paste("Lambda =",lambda,"\n"))
	##working from metal output###
	if (is.null(sig)==TRUE){
		ymin1=round( max(-log10(data$P.Value))+1)
	} else {
		ymin1=12
	}
	title=c()
	jpeg(paste("manhattan_", filename, ".jpeg",sep=""),res=400,width = 40, height = 12,units="cm")
	chr <- c(1:22)
	#Summary statistics
	data$position<-round(data$start,digits=0)
	#print(summary(data))
	print(table(data$chr))
	par(mar=c(5,5,2,2))
	phy.max<-tapply(data$start, data$chr,max,na.rm=T)
	cumlen=0
	data <- data[!is.na(data$chr),]
	for(i in chr){
		cat(paste("Now working on chromosome ",i,"\r"))
		data[data$chr==i,"loc"]<-data[data$chr==i,"position"]+cumlen
		cumlen<-cumlen+phy.max[i]
	}
	phy.med<-tapply(data$loc,data$chr,median,na.rm=T)[chr]
	print(phy.med)
	data$mlgpval<- -log(data[,"P.Value"], base=10)
	plot(data[,"loc"],data[,"mlgpval"],type="n",yaxt="n",xaxt="n",
		xlab="chromosome",
		ylab=expression(-log[10]*P),main=title,
		xlim=c(0,max(data$loc,na.rm=T)),cex.lab=1.5,ylim=c(0,ymin1))
		col=ifelse(data$mlgpval< 10, (rep(c("black","gray48"),13)), "red" )
	axis(side=2, at=seq(from=0,to=ymin1,by=2), labels=seq(from=0,to=ymin1,by=2),
		tick=T,cex.axis=0.9,las=1)
	axis(side=1, at=phy.med[c(1:22)], labels=chr[c(1:22)],
		tick=T,cex.axis=0.7,las=1)

		#cex.axis=0.8 instead of cex.axis=1.2
	#axis(side=1, at=phy.med[c(20:23)], labels=chr[c(20:23)],
	#	tick=T,cex.axis=1,las=1)
	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
		col = "white")
	for(i in chr){
		cat(paste("Now working on chromosome ",i,"\r"))
		#if(data$mlgpval > 7) col="blue" else col="red"
		points(data[data$chr==i,"loc"],data[data$chr==i,"mlgpval"],
			col=ifelse(data[data$chr==i,"mlgpval"] > 10, "red", col[i]),pch=20) #,cex=data[data$chr==i,"cex.val"])
		#,cex = data[data$chr==i,"cex.val"]
	}
	# add in line for significance
	if (is.null(sig)==FALSE){
		abline(h=sig,lty="dotted",lwd=2,col="chartreuse4")
	}

	dev.off()
}
################

manhattan(DMS=DMS, filename="CCForPublicationEPIC", sig=7.3)

###############################################################

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

# code to find differentially methylated regions

library(DMRcate)

type <- factor(covariates$Line)
design <- model.matrix(~type+Batch+ID)
myannotation <- cpg.annotate("array", EPICTFAMKO, arraytype="EPIC", analysis.type="differential", design=design, coef=2, what="Beta", fdr=0.01)

dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2, min.cpgs=10)

results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
results <- as.data.frame(results.ranges)

groups <- c(KO="magenta", NC="forestgreen")
cols <- groups[as.character(type)]