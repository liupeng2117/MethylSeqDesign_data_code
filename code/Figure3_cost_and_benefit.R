library(pi0) #estimate pi0(lambda)
library(ggplot2) 
library(VGAM) # estimate fai
library(fitdistrplus) # fit empirical distribution
library(doParallel) # foreach
library(knitr)
library(plyr) #map values of pilotN based on key

setwd("D:/research/MethySeq/")
source("code/functions.R")
load(file="data/Mouse/meth.regional.count.R")
load(file="data/Mouse/empirical.fit.rdata")

pilot.N=4
N<-seq(4,50,by=2)
pilot.R=1/4
R=c(1/10, 1/8, 1/6, 1/4)
#prop=(R/pilot.R)[(R/pilot.R)<=1]
rep<-10
G<-10000
n.DE<-1000
set.seed(123)
delta<-runif(n.DE,0.1,0.2)
set.seed(123)
dmr<-sample(1:G,n.DE)


simulated.result<-simulate.methyl.data(n=pilot.N, G=G, delta=delta, dmr=dmr, prop=1)
dmr.result<-DMR.analysis(N0=pilot.N,cov.matrix=simulated.result$coverage, methyl.matrix=simulated.result$methyl.count, R=R, pilot.depth=250/8)
EDR.result<-Estimate.EDR.from.pilot(res=dmr.result, thresh.p=0.005,
                                       N0=pilot.N, target.N=N, FDR=0.05, M=10)


plot2d(EDR.result)
plot3d(EDR.result)

design.case1<-designOptim(EDR.result, pilot.depth=250/8, R, N, targetEDR=0.8)
design.case2<-designOptim(EDR.result, pilot.depth=250/8, R, N, budget=20000)

library(gridExtra)
pdf("results/cost and benefit.pdf")
grid.arrange(design.case1$plot1,design.case2$plot1, design.case1$plot2, design.case2$plot2, ncol=2)
dev.off()

