library(ggplot2) 
library(VGAM) # estimate fai
library(fitdistrplus) # fit empirical distribution
library(methylKit) # read in the raw data and do region based analysis
library(doParallel) # foreach
library(knitr)
library(plyr) #map values of pilotN based on key



## 1. Region based methylation dataset
setwd("D:/research/MethySeq/")
source("code/functions.R")
load(file="data/Mouse/meth.regional.count.R")

## 2.Simulate region-based methylation dataset

##########################################
# Here we need to estimate a common dispersion first
# Then estimate empirical distribution of coverage,draw mean from the distribution
# Then estimate size for negative binomial distribution
# Using mean drawn and estimated size, randomly generate coverage
# Then estimate empirical distribution of p, draw mean of p
# With common dispersion and mean of p, calculate alpha1 and beta1
# Simulate the signal:
  # Set a 0<delta<0.5
  # If p>0.5, p'=p-delta
  # If p<0.5, p'=p+delta
# With common dispersion and mean of p', calculate alpha2 and beta2 
# With alpha1, beta1 and alpha2,beta2, simulate p for each combination of region and sample
# With p and coverage, use binomial model to simulate number of methylated count
##########################################

### 2.1 parameters setting
#### 2.1.1 Here first estimate common dispersion parameter fai for beta distribution of p
    # To save time I ran fist 1000
seq.coverage<-seq(5,32,by=3)
seq.numC<-seq(6,33,by=3)
fai.mean<-0.04838581

#### 2.1.2 Then estimate the mean dispersion parameter(size) for negative binomial distribuiton of coverage
coverage<-meth[,seq.coverage]
# estimate the dispersion
#size.mean<-mean(size) #6.286579
 size.mean<-6.286579

#### 2.1.3 Fit empirical distribution of coverage
coverage<-meth[,seq.coverage]
coverage.mean<-apply(coverage,1,mean)
fitemp.c<-ecdf(coverage.mean)


#### 2.1.4 Fit empirical distributio of p

    # Estimate p from real data
p.estimate<-matrix(0,nrow(meth),5)
for (i in 6:10)
{
  numC.col<-(6+3*(i-1))
  coverage.col<-(5+3*(i-1))
  p.estimate[,(i-5)]<-meth[,numC.col]/meth[,coverage.col]
}

p.estimate.mean<-apply(p.estimate,1,mean)
fitemp.p<-ecdf(p.estimate.mean)


#### 2.1.5 Number of samples, regions and DMRs in the simulated dataset
D<-10
G<-10000
n.DE<-1000
set.seed(123)
dmr<-sample(1:G,n.DE)


## 3. Comparing 4 tests
#####################
# Beta value + t test
# M value + t test
# Z value + t test
# Z value + wald test
# beta value= #methy/#coverage
# M value=logit(beta)
# Z value=arcsin(2*beta-1)
#####################
iter<-20
tpr.n1<-matrix(0,iter,n.DE)
tpr.r1<-matrix(0,iter,n.DE)
tpr.n2<-matrix(0,iter,n.DE)
tpr.r2<-matrix(0,iter,n.DE)
tpr.n3<-matrix(0,iter,n.DE)
tpr.r3<-matrix(0,iter,n.DE)
tpr.r4<-matrix(0,iter,n.DE)
tpr.n4<-matrix(0,iter,n.DE)

tpr.n1.low<-matrix(0,iter,n.DE)
tpr.r1.low<-matrix(0,iter,n.DE)
tpr.n2.low<-matrix(0,iter,n.DE)
tpr.r2.low<-matrix(0,iter,n.DE)
tpr.n3.low<-matrix(0,iter,n.DE)
tpr.r3.low<-matrix(0,iter,n.DE)
tpr.n4.low<-matrix(0,iter,n.DE)
tpr.r4.low<-matrix(0,iter,n.DE)

tpr.n1.medium<-matrix(0,iter,n.DE)
tpr.r1.medium<-matrix(0,iter,n.DE)
tpr.n2.medium<-matrix(0,iter,n.DE)
tpr.r2.medium<-matrix(0,iter,n.DE)
tpr.n3.medium<-matrix(0,iter,n.DE)
tpr.r3.medium<-matrix(0,iter,n.DE)
tpr.n4.medium<-matrix(0,iter,n.DE)
tpr.r4.medium<-matrix(0,iter,n.DE)

tpr.n1.high<-matrix(0,iter,n.DE)
tpr.r1.high<-matrix(0,iter,n.DE)
tpr.n2.high<-matrix(0,iter,n.DE)
tpr.r2.high<-matrix(0,iter,n.DE)
tpr.n3.high<-matrix(0,iter,n.DE)
tpr.r3.high<-matrix(0,iter,n.DE)
tpr.n4.high<-matrix(0,iter,n.DE)
tpr.r4.high<-matrix(0,iter,n.DE)

### 3.1 Overall 
for (t in 1:iter){ # 20 iterations
  simulated.result<-simulate.methyl.data(D=D,G=G,delta=0.15,dmr=dmr)
  # p 
  p.Bt<-B.ttest(sim.coverage=simulated.result$coverage,sim.meth=simulated.result$methyl.count,D=D)
  p.Mt<-M.ttest(sim.coverage=simulated.result$coverage,sim.meth=simulated.result$methyl.count,D=D)
  p.Zt<-Z.ttest(sim.coverage=simulated.result$coverage,sim.meth=simulated.result$methyl.count,D=D)
  fai.est<-estimate.fai2(sim.coverage = simulated.result$coverage,sim.meth = simulated.result$methyl.count)
  test.result<-Z.wald(sim.coverage=simulated.result$coverage,sim.meth=simulated.result$methyl.count,fai.est=fai.est,D=D)
  p.Zw<-test.result[[1]]
  # p.value Beta+t
  p.Bt<-cbind(region=1:length(p.Bt),p=p.Bt)
  p.Bt<-p.Bt[order(p.Bt[,"p"]),]
  # ture positives
  for(i in 1:n.DE){ # top n.DE p values 
    tpr.n1[t,i]<-length(intersect(p.Bt[1:i,1],dmr))
  }  
  tpr.r1[t,]<-tpr.n1[t,]/1:n.DE
  
  
  # p.value M+t
  p.Mt<-cbind(region=1:length(p.Mt),p=p.Mt)
  p.Mt<-p.Mt[order(p.Mt[,"p"]),]
  # ture positives
  for(i in 1:n.DE){
    tpr.n2[t,i]<-length(intersect(p.Mt[1:i,1],dmr))
  }  
  tpr.r2[t,]<-tpr.n2[t,]/1:n.DE
  
  # p.value Z+t
  p.Zt<-cbind(region=1:length(p.Zt),p=p.Zt)
  p.Zt<-p.Zt[order(p.Zt[,"p"]),]
  #true positives
  for(i in 1:n.DE){
    tpr.n3[t,i]<-length(intersect(p.Zt[1:i,1],dmr))
  }  
  tpr.r3[t,]<-tpr.n3[t,]/1:n.DE
  
  # p.value Z+w
  p.Zw<-cbind(region=1:length(p.Zw),p=p.Zw)
  p.Zw<-p.Zw[order(p.Zw[,"p"]),]
  #true positives
  for(i in 1:n.DE){
    tpr.n4[t,i]<-length(intersect(p.Zw[1:i,1],dmr))
  }  
  tpr.r4[t,]<-tpr.n4[t,]/1:n.DE
  }

### 3.2 Low
for (t in 1:iter){
  simulated.result<-simulate.methyl.data(D=D,G=G,delta=0.15,up.bound = 0.2,dmr=dmr)
  # p
  p.Bt.low<-B.ttest(sim.coverage=simulated.result$coverage,sim.meth =simulated.result$methyl.count,D=D)
  p.Mt.low<-M.ttest(sim.coverage=simulated.result$coverage,sim.meth=simulated.result$methyl.count,D=D)
  p.Zt.low<-Z.ttest(sim.coverage=simulated.result$coverage,sim.meth=simulated.result$methyl.count,D=D)
  fai.est<-estimate.fai2(sim.coverage = simulated.result$coverage,sim.meth = simulated.result$methyl.count)
  test.result.low<-Z.wald(sim.coverage=simulated.result$coverage,sim.meth=simulated.result$methyl.count,fai.est=fai.est,D=D)
  p.Zw.low<-test.result.low[[1]]
  # p.value beta+t low
  p.Bt.low<-cbind(region=1:length(p.Bt.low),p=p.Bt.low)
  p.Bt.low<-p.Bt.low[order(p.Bt.low[,"p"]),]
  # true positive
  for(i in 1:n.DE){
    tpr.n1.low[t,i]<-length(intersect(p.Bt.low[1:i,1],dmr))
  }  
  tpr.r1.low[t,]<-tpr.n1.low[t,]/1:n.DE
  
  # p.value M+t low
  p.Mt.low<-cbind(region=1:length(p.Mt.low),p=p.Mt.low)
  p.Mt.low<-p.Mt.low[order(p.Mt.low[,"p"]),]
  # true positive
  for(i in 1:n.DE){
    tpr.n2.low[t,i]<-length(intersect(p.Mt.low[1:i,1],dmr))
  }  
  tpr.r2.low[t,]<-tpr.n2.low[t,]/1:n.DE
  
  # p.value Z+t low
  p.Zt.low<-cbind(region=1:length(p.Zt.low),p=p.Zt.low)
  p.Zt.low<-p.Zt.low[order(p.Zt.low[,"p"]),]
  # true positive
  for(i in 1:n.DE){
    tpr.n3.low[t,i]<-length(intersect(p.Zt.low[1:i,1],dmr))
  }  
  tpr.r3.low[t,]<-tpr.n3.low[t,]/1:n.DE
  
  # p.value Z+w low
  p.Zw.low<-cbind(region=1:length(p.Zw.low),p=p.Zw.low)
  p.Zw.low<-p.Zw.low[order(p.Zw.low[,"p"]),]
  # true positive
  for(i in 1:n.DE){
    tpr.n4.low[t,i]<-length(intersect(p.Zw.low[1:i,1],dmr))
  }  
  tpr.r4.low[t,]<-tpr.n4.low[t,]/1:n.DE
}

### 3.3 medium
for (t in 1:iter){
  simulated.result<-simulate.methyl.data(D=D,G=G,delta=0.15,low.bound = 0.2,up.bound = 0.8,dmr=dmr)
  
  # p values
  p.Bt.medium<-B.ttest(sim.coverage=simulated.result$coverage,sim.meth=simulated.result$methyl.count,D=D)
  p.Mt.medium<-M.ttest(sim.coverage=simulated.result$coverage,sim.meth=simulated.result$methyl.count,D=D)
  p.Zt.medium<-Z.ttest(sim.coverage=simulated.result$coverage,sim.meth=simulated.result$methyl.count,D=D)
  fai.est<-estimate.fai2(sim.coverage = simulated.result$coverage,sim.meth = simulated.result$methyl.count)
  test.result.medium<-Z.wald(sim.coverage=simulated.result$coverage,sim.meth=simulated.result$methyl.count,fai.est=fai.est,D=D)
  p.Zw.medium<-test.result.medium[[1]]
  # p.value B+t medium
  p.Bt.medium<-cbind(region=1:length(p.Bt.medium),p=p.Bt.medium)
  p.Bt.medium<-p.Bt.medium[order(p.Bt.medium[,"p"]),]
  # true positive
  for(i in 1:n.DE){
    tpr.n1.medium[t,i]<-length(intersect(p.Bt.medium[1:i,1],dmr))
  }  
  tpr.r1.medium[t,]<-tpr.n1.medium[t,]/1:n.DE
  
  # p.value M+t medium
  p.Mt.medium<-cbind(region=1:length(p.Mt.medium),p=p.Mt.medium)
  p.Mt.medium<-p.Mt.medium[order(p.Mt.medium[,"p"]),]
  # true positive
  for(i in 1:n.DE){
    tpr.n2.medium[t,i]<-length(intersect(p.Mt.medium[1:i,1],dmr))
  }  
  tpr.r2.medium[t,]<-tpr.n2.medium[t,]/1:n.DE
  
  # p.value Z+t medium
  p.Zt.medium<-cbind(region=1:length(p.Zt.medium),p=p.Zt.medium)
  p.Zt.medium<-p.Zt.medium[order(p.Zt.medium[,"p"]),]
  # true positive
  for(i in 1:n.DE){
    tpr.n3.medium[t,i]<-length(intersect(p.Zt.medium[1:i,1],dmr))
  }  
  tpr.r3.medium[t,]<-tpr.n3.medium[t,]/1:n.DE
  
  # p.value Z+w medium
  p.Zw.medium<-cbind(region=1:length(p.Zw.medium),p=p.Zw.medium)
  p.Zw.medium<-p.Zw.medium[order(p.Zw.medium[,"p"]),]
  # true positive
  for(i in 1:n.DE){
    tpr.n4.medium[t,i]<-length(intersect(p.Zw.medium[1:i,1],dmr))
  }  
  tpr.r4.medium[t,]<-tpr.n4.medium[t,]/1:n.DE
}

### 3.4 High
for (t in 1:iter){
  simulated.result<-simulate.methyl.data(D=D,G=G,delta=0.15,low.bound = 0.8,dmr=dmr)
  
  # p values
  p.Bt.high<-B.ttest(sim.coverage=simulated.result$coverage,sim.meth=simulated.result$methyl.count,D=D)
  p.Mt.high<-M.ttest(sim.coverage=simulated.result$coverage,sim.meth=simulated.result$methyl.count,D=D)
  p.Zt.high<-Z.ttest(sim.coverage=simulated.result$coverage,sim.meth=simulated.result$methyl.count,D=D)
  fai.est<-estimate.fai2(sim.coverage = simulated.result$coverage,sim.meth = simulated.result$methyl.count)
  test.result<-Z.wald(sim.coverage=simulated.result$coverage,sim.meth=simulated.result$methyl.count,fai.est=fai.est,D=D)
  p.Zw.high<-test.result[[1]]
  # p.value M+t high
  p.Bt.high<-cbind(region=1:length(p.Bt.high),p=p.Bt.high)
  p.Bt.high<-p.Bt.high[order(p.Bt.high[,"p"]),]
  # true positive
  for(i in 1:n.DE){
    tpr.n1.high[t,i]<-length(intersect(p.Bt.high[1:i,1],dmr))
  }  
  tpr.r1.high[t,]<-tpr.n1.high[t,]/1:n.DE
  
  # p.value M+t high
  p.Mt.high<-cbind(region=1:length(p.Mt.high),p=p.Mt.high)
  p.Mt.high<-p.Mt.high[order(p.Mt.high[,"p"]),]
  # true positive
  for(i in 1:n.DE){
    tpr.n2.high[t,i]<-length(intersect(p.Mt.high[1:i,1],dmr))
  }  
  tpr.r2.high[t,]<-tpr.n2.high[t,]/1:n.DE
  
  # p.value Z+t high
  p.Zt.high<-cbind(region=1:length(p.Zt.high),p=p.Zt.high)
  p.Zt.high<-p.Zt.high[order(p.Zt.high[,"p"]),]
  # true positive
  for(i in 1:n.DE){
    tpr.n3.high[t,i]<-length(intersect(p.Zt.high[1:i,1],dmr))
  }  
  tpr.r3.high[t,]<-tpr.n3.high[t,]/1:n.DE
  
  # p.value Z+w high
  p.Zw.high<-cbind(region=1:length(p.Zw.high),p=p.Zw.high)
  p.Zw.high<-p.Zw.high[order(p.Zw.high[,"p"]),]
  # true positive
  for(i in 1:n.DE){
    tpr.n4.high[t,i]<-length(intersect(p.Zw.high[1:i,1],dmr))
  }  
  tpr.r4.high[t,]<-tpr.n4.high[t,]/1:n.DE
}

save.image("data/simulation/Simulation_stratification.Rdata")
##########plot#################
  #plot1
  pdf("results/straitification2.pdf", width = 8, height = 8)
  par(mfrow=c(2,2))
  plot(1:n.DE,apply(tpr.n1,2,mean),type="l",lty=3,lwd=2,col=8,ylim=c(1,n.DE),
       xlab="The number of most significant regions selected",ylab="The number of true positives",
       main="Overall")
  lines(c(1,n.DE),c(1,n.DE),col="red")
  par(new=T)
  plot(1:n.DE,apply(tpr.n3,2,mean),type="l",lwd=2,ylim=c(1,n.DE),
       xlab="",ylab="",axes=F,col=8)
  
  par(new=T)
  plot(1:n.DE,apply(tpr.n2,2,mean),type="l",lty=5,lwd=2,ylim=c(1,n.DE),
       xlab="",ylab="",axes=F)
  
   par(new=T)
  plot(1:n.DE,apply(tpr.n4,2,mean),type="l",lwd=2,ylim=c(1,n.DE),
       xlab="",ylab="",axes=F)
  
    # low plot1
plot(1:n.DE,apply(tpr.n1.low,2,mean),type="l",lty=3,lwd=2,col=8,ylim=c(1,n.DE),
     xlab="The number of most significant regions selected",ylab="The number of true positives",
     main="Low")
lines(c(1,n.DE),c(1,n.DE),col="red")
par(new=T)
plot(1:n.DE,apply(tpr.n3.low,2,mean),type="l",lwd=2,ylim=c(1,n.DE),
     xlab="",ylab="",axes=F,col=8)

par(new=T)
plot(1:n.DE,apply(tpr.n2.low,2,mean),type="l",lty=5,lwd=2,ylim=c(1,n.DE),
     xlab="",ylab="",axes=F)
par(new=T)
plot(1:n.DE,apply(tpr.n4.low,2,mean),type="l",lwd=2,ylim=c(1,n.DE),
     xlab="",ylab="",axes=F)

   #medium plot1 
plot(1:n.DE,apply(tpr.n1.medium,2,mean),type="l",lty=3,lwd=2,col=8,ylim=c(1,n.DE),
     xlab="The number of most significant regions selected",ylab="The number of true positives",
     main="medium")
lines(c(1,n.DE),c(1,n.DE),col="red")
par(new=T)
plot(1:n.DE,apply(tpr.n3.medium,2,mean),type="l",lwd=2,ylim=c(1,n.DE),
     xlab="",ylab="",axes=F,col=8)

par(new=T)
plot(1:n.DE,apply(tpr.n2.medium,2,mean),type="l",lty=5,lwd=2,ylim=c(1,n.DE),
     xlab="",ylab="",axes=F)
par(new=T)
plot(1:n.DE,apply(tpr.n4.medium,2,mean),type="l",lwd=2,ylim=c(1,n.DE),
     xlab="",ylab="",axes=F)

    #high plot1 
plot(1:n.DE,apply(tpr.n1.high,2,mean),type="l",lty=3,lwd=2,col=8,ylim=c(1,n.DE),
     xlab="The number of most significant regions selected",ylab="The number of true positives",
     main="high")
lines(c(1,n.DE),c(1,n.DE),col="red")
par(new=T)
plot(1:n.DE,apply(tpr.n3.high,2,mean),type="l",lwd=2,ylim=c(1,n.DE),
     xlab="",ylab="",axes=F,col=8)

par(new=T)
plot(1:n.DE,apply(tpr.n2.high,2,mean),type="l",lty=5,lwd=2,ylim=c(1,n.DE),
     xlab="",ylab="",axes=F)
par(new=T)
plot(1:n.DE,apply(tpr.n4.high,2,mean),type="l",lwd=2,ylim=c(1,n.DE),
     xlab="",ylab="",axes=F)
legend( x="bottomright", 
        legend=c("A+Wald","A+t","M+t","Beta+t"),
        col=c(1,8,1,8), lwd=2, lty=c(1,1,5,3))

dev.off()
