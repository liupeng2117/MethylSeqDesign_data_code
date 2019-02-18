setwd("D:/research/MethySeq/")
source("functions.R")
library(ggplot2) 
library(VGAM) # estimate fai
library(fitdistrplus) # fit empirical distribution
library(methylKit) # read in the raw data and do region based analysis
library(doParallel) # foreach
library(knitr)
library(plyr) #map values of pilotN based on key

## 1. Region based methylation dataset
setwd("D:/research/MethySeq/Mouse")
bead.info <- read.table("WUSTL_BaitTiling_merged.bed")
ROI.regions = cbind(bead.info,rep(NA,nrow(bead.info)))
colnames(ROI.regions) = c("chr","start","end","seqnames")
require(GenomicRanges)

data.frame2GRanges <- function(df, keepColumns = FALSE, ignoreStrand = FALSE) {
  stopifnot(class(df) == "data.frame")
  stopifnot(all(c("start", "end") %in% names(df)))
  stopifnot(any(c("chr", "seqnames") %in% names(df)))
  if("seqnames" %in% names(df))
    names(df)[names(df) == "seqnames"] <- "chr"
  if(!ignoreStrand && "strand" %in% names(df)) {
    if(is.numeric(df$strand)) {
      strand <- ifelse(df$strand == 1, "+", "*")
      strand[df$strand == -1] <- "-"
      df$strand <- strand
    }
    gr <- GRanges(seqnames = df$chr,
                  ranges = IRanges(start = df$start, end = df$end),
                  strand = df$strand)
  } else {
    gr <- GRanges(seqnames = df$chr,
                  ranges = IRanges(start = df$start, end = df$end))
  }
  if(keepColumns) {
    dt <- as(df[, setdiff(names(df), c("chr", "start", "end", "strand"))],
             "DataFrame")
    elementMetadata(gr) <- dt
  }
  names(gr) <- rownames(df)
  gr
}

ROI.Granges = data.frame2GRanges(ROI.regions)


my.path<-"D:/research/MethySeq/Mouse/MethylKit_P_1_batch_1"  # the data is publicly available on https://figshare.com/articles/MethylKit_P_1_batch_1/7730753
setwd(my.path)
file.names = sort(list.files(pattern="_CpG.txt"))
classlabel = rep(0,length(file.names))
classlabel[which(substr(file.names,1,1)%in%"V")] = 0
classlabel[which(substr(file.names,1,1)%in%"P")] = 1


sample.id = lapply(strsplit(file.names,"_"),function(x) x[1])
file.names.new = as.list(file.names)
classlabel.new = classlabel[-1]
myobj <- methRead(file.names.new, sample.id = sample.id, assembly = "mm9", treatment = classlabel.new, context = "CpG")

# Get the count per region
regional.count<-regionCounts(myobj,regions=ROI.Granges)

## filtering
# filter out counts<20 and >99.9th percentile
filtered.regional.count=filterByCoverage(regional.count,lo.count=20,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)

#  merge all samples to one object for regions covered in all samples
meth.regional.count <- unite(filtered.regional.count)
dim(meth.regional.count) ## 137148 regions

meth<-getData(meth.regional.count)
save(meth,file="meth.regional.count.R")

#estimate fai
seq.coverage<-seq(5,32,by=3)
seq.numC<-seq(6,33,by=3)
coverage<-meth[,seq.coverage]
numC<-meth[,seq.numC]
coverage<-as.matrix(coverage)
numC<-as.matrix(numC)
#Use VGAM
fai<-sapply(1:nrow(numc),function(x) {fit=vglm(cbind(numC[x,],coverage[x,]-numC[x,])~1, betabinomial)
                           c(mean(fit@misc$rho),mean(fit@misc$size))})
fai.mean=mean(fai)
#Use estimate.fai
fai.mean<-estimate.fai(coverage, numC, x=c(5,5))
save(fai.mean,file="empirical.fit.rdata")