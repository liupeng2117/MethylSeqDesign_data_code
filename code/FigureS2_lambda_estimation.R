library(pi0) #estimate pi0(lambda)
library(ggplot2) 
library(gridExtra)
library(grid)
library(VGAM) # estimate fai
library(fitdistrplus) # fit empirical distribution
library(doParallel) # foreach
library(knitr)
library(plyr) #map values of pilotN based on key

setwd("D:/research/MethySeq/")
source("code/functions.R")
load(file="data/Mouse/meth.regional.count.R")
load(file="data/Mouse/empirical.fit.rdata")

pilot.N.all<-c(4, 8, 12, 16, 18, 20)
N<-c(2,6,10,15,25,50)
rep<-10
G<-10000
n.DE<-1000
set.seed(123)
dmr<-sample(1:G,n.DE)

#--------------------------
# delta=0.12
#--------------------------
delta<-rep(0.12,n.DE)

#Assign 8 cores
cl<-makeCluster(8)
registerDoParallel(cl)
result<-vector("list",length(pilot.N.all))

lambda.allN0<-list()
i=1
for(pilot.N in pilot.N.all){
  lambda.oneN0<-foreach (times = 1:rep,.packages = "pi0") %dopar% {
    #simulate data
    simulated.result<-simulate.methyl.data(n=pilot.N, G=G, delta=delta, dmr=dmr)
    
    #DE analysis
    a=proc.time()
    dmr.result<-DMR.analysis(N0=pilot.N,cov.matrix=simulated.result$coverage, methyl.matrix=simulated.result$methyl.count, R=1/4, pilot.R=1/4)
    #estimate pi0
    estimated.lambda<-Estimate.lambda.from.pilot(p.values=dmr.result$p.values, thresh.p=0.005,
                                                 N0 = pilot.N, target.N = N, FDR=0.05, M=10)
    b=proc.time()-a
    return(c(estimated.lambda$lambda, estimated.lambda$lambda.CDD, estimated.lambda$lambda.MLE))
  }
  lambda.oneN0<-Reduce(rbind,lambda.oneN0)
  colnames(lambda.oneN0)<- c("CBUM", "CDD", "MLE")
  lambda.allN0[[i]]<-lambda.oneN0
  i=i+1
}

stopCluster(cl)

names(lambda.allN0)<-as.character(pilot.N.all)
save(lambda.allN0,file="data/simulation/lambda.estimation.delta0.12.cleaned.rdata")

lambda.mean<-sapply(lambda.allN0, function(x) mean(x[,1]))
lambda.l<-sapply(lambda.allN0, function(x) quantile(x[,1],probs=0.25))
lambda.u<-sapply(lambda.allN0, function(x) quantile(x[,1],probs=0.75))
cbum.lambda<-data.frame(method="cbum",n0=pilot.N.all/2, mean=lambda.mean, low=lambda.l, up=lambda.u)

lambda.mean<-sapply(lambda.allN0, function(x) mean(x[,3]))
lambda.l<-sapply(lambda.allN0, function(x) quantile(x[,3],probs=0.25))
lambda.u<-sapply(lambda.allN0, function(x) quantile(x[,3],probs=0.75))
mle.lambda<-data.frame(method="mle",n0=pilot.N.all/2, mean=lambda.mean, low=lambda.l, up=lambda.u)

lambda.data.0.12<-rbind(cbum.lambda, mle.lambda)
#--------------------------
# delta=0.15
#--------------------------

delta<-rep(0.15,n.DE)

#Assign 8 cores
cl<-makeCluster(8)
registerDoParallel(cl)
result<-vector("list",length(pilot.N.all))

lambda.allN0<-list()
i=1
for(pilot.N in pilot.N.all){
  lambda.oneN0<-foreach (times = 1:rep,.packages = "pi0") %dopar% {
    #simulate data
    simulated.result<-simulate.methyl.data(n=pilot.N, G=G, delta=delta, dmr=dmr)
    
    #DE analysis
    a=proc.time()
    dmr.result<-DMR.analysis(N0=pilot.N,cov.matrix=simulated.result$coverage, methyl.matrix=simulated.result$methyl.count, R=1/4, pilot.R=1/4)
    #estimate pi0
    estimated.lambda<-Estimate.lambda.from.pilot(p.values=dmr.result$p.values, thresh.p=0.005,
                                                 N0 = pilot.N, target.N = N, FDR=0.05, M=10)
    b=proc.time()-a
    return(c(estimated.lambda$lambda, estimated.lambda$lambda.CDD, estimated.lambda$lambda.MLE))
  }
  lambda.oneN0<-Reduce(rbind,lambda.oneN0)
  colnames(lambda.oneN0)<- c("CBUM", "CDD", "MLE")
  lambda.allN0[[i]]<-lambda.oneN0
  i=i+1
}
stopCluster(cl)

names(lambda.allN0)<-as.character(pilot.N.all)
save(lambda.allN0,file="data/simulation/lambda.estimation.delta0.15.cleaned.rdata")

lambda.mean<-sapply(lambda.allN0, function(x) mean(x[,1]))
lambda.l<-sapply(lambda.allN0, function(x) quantile(x[,1],probs=0.25))
lambda.u<-sapply(lambda.allN0, function(x) quantile(x[,1],probs=0.75))
cbum.lambda<-data.frame(method="cbum",n0=pilot.N.all/2, mean=lambda.mean, low=lambda.l, up=lambda.u)

lambda.mean<-sapply(lambda.allN0, function(x) mean(x[,3]))
lambda.l<-sapply(lambda.allN0, function(x) quantile(x[,3],probs=0.25))
lambda.u<-sapply(lambda.allN0, function(x) quantile(x[,3],probs=0.75))
mle.lambda<-data.frame(method="mle",n0=pilot.N.all/2, mean=lambda.mean, low=lambda.l, up=lambda.u)

lambda.data.0.15<-rbind(cbum.lambda, mle.lambda)

#--------------------------
# delta=0.18
#--------------------------

delta<-rep(0.18,n.DE)

#Assign 8 cores
cl<-makeCluster(8)
registerDoParallel(cl)
result<-vector("list",length(pilot.N.all))

lambda.allN0<-list()
i=1
for(pilot.N in pilot.N.all){
  lambda.oneN0<-foreach (times = 1:rep,.packages = "pi0") %dopar% {
    #simulate data
    simulated.result<-simulate.methyl.data(n=pilot.N, G=G, delta=delta, dmr=dmr)
    
    #DE analysis
    a=proc.time()
    dmr.result<-DMR.analysis(N0=pilot.N,cov.matrix=simulated.result$coverage, methyl.matrix=simulated.result$methyl.count, R=1/4, pilot.R=1/4)
    #estimate pi0
    estimated.lambda<-Estimate.lambda.from.pilot(p.values=dmr.result$p.values, thresh.p=0.005,
                                                 N0 = pilot.N, target.N = N, FDR=0.05, M=10)
    b=proc.time()-a
    return(c(estimated.lambda$lambda, estimated.lambda$lambda.CDD, estimated.lambda$lambda.MLE))
  }
  lambda.oneN0<-Reduce(rbind,lambda.oneN0)
  colnames(lambda.oneN0)<- c("CBUM", "CDD", "MLE")
  lambda.allN0[[i]]<-lambda.oneN0
  i=i+1
}
stopCluster(cl)

names(lambda.allN0)<-as.character(pilot.N.all)
save(lambda.allN0,file="data/simulation/lambda.estimation.delta0.18.cleaned.rdata")

lambda.mean<-sapply(lambda.allN0, function(x) mean(x[,1]))
lambda.l<-sapply(lambda.allN0, function(x) quantile(x[,1],probs=0.25))
lambda.u<-sapply(lambda.allN0, function(x) quantile(x[,1],probs=0.75))
cbum.lambda<-data.frame(method="cbum",n0=pilot.N.all/2, mean=lambda.mean, low=lambda.l, up=lambda.u)

lambda.mean<-sapply(lambda.allN0, function(x) mean(x[,3]))
lambda.l<-sapply(lambda.allN0, function(x) quantile(x[,3],probs=0.25))
lambda.u<-sapply(lambda.allN0, function(x) quantile(x[,3],probs=0.75))
mle.lambda<-data.frame(method="mle",n0=pilot.N.all/2, mean=lambda.mean, low=lambda.l, up=lambda.u)

lambda.data.0.18<-rbind(cbum.lambda, mle.lambda)


#--------------------------
# uniform delta
#--------------------------
delta<-runif(n.DE,0.1,0.2)

#Assign 8 cores
cl<-makeCluster(8)
registerDoParallel(cl)
result<-vector("list",length(pilot.N.all))

lambda.allN0<-list()
i=1
for(pilot.N in pilot.N.all){
      lambda.oneN0<-foreach (times = 1:rep,.packages = "pi0") %dopar% {
        #simulate data
        simulated.result<-simulate.methyl.data(n=pilot.N, G=G, delta=delta, dmr=dmr)

        #DE analysis
        a=proc.time()
        dmr.result<-DMR.analysis(N0=pilot.N,cov.matrix=simulated.result$coverage, methyl.matrix=simulated.result$methyl.count,R=1/4, pilot.R=1/4)
        #estimate pi0
        estimated.lambda<-Estimate.lambda.from.pilot(p.values=dmr.result$p.values, thresh.p=0.005,
                                              N0 = pilot.N, target.N = N, FDR=0.05, M=10)
        b=proc.time()-a
        return(c(estimated.lambda$lambda, estimated.lambda$lambda.CDD, estimated.lambda$lambda.MLE))
      }
    lambda.oneN0<-Reduce(rbind,lambda.oneN0)
    colnames(lambda.oneN0)<- c("CBUM", "CDD", "MLE")
    lambda.allN0[[i]]<-lambda.oneN0
    i=i+1
}
stopCluster(cl)

names(lambda.allN0)<-as.character(pilot.N.all)
save(lambda.allN0,file="data/simulation/lambda.estimation.unif.cleaned.rdata")

lambda.mean<-sapply(lambda.allN0, function(x) mean(x[,1]))
lambda.l<-sapply(lambda.allN0, function(x) quantile(x[,1],probs=0.25))
lambda.u<-sapply(lambda.allN0, function(x) quantile(x[,1],probs=0.75))
cbum.lambda<-data.frame(method="cbum",n0=pilot.N.all/2, mean=lambda.mean, low=lambda.l, up=lambda.u)

lambda.mean<-sapply(lambda.allN0, function(x) mean(x[,3]))
lambda.l<-sapply(lambda.allN0, function(x) quantile(x[,3],probs=0.25))
lambda.u<-sapply(lambda.allN0, function(x) quantile(x[,3],probs=0.75))
mle.lambda<-data.frame(method="mle",n0=pilot.N.all/2, mean=lambda.mean, low=lambda.l, up=lambda.u)

lambda.data.unif<-rbind(cbum.lambda, mle.lambda)

#------------------
# plot
#------------------

p1<-ggplot(lambda.data.0.12, aes(x=n0, y=mean)) + 
  geom_errorbar(aes(ymin=low, ymax=up), width=0.2) +
  geom_line(lwd=1, aes(linetype=method)) +
  geom_point()+
  ylim(0.7,1)+xlim(1,11)+
  labs(x = "Sample size per group", y = "Estimated lambda", title="Delta = 0.12") +
  scale_x_continuous(breaks=seq(0,12,2))+
  geom_hline(yintercept=0.9, linetype="dashed", color = "red")+
  theme_bw()+theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

p2<-ggplot(lambda.data.0.15, aes(x=n0, y=mean)) + 
  geom_errorbar(aes(ymin=low, ymax=up), width=0.2) +
  geom_line(lwd=1, aes(linetype=method)) +
  geom_point()+
  ylim(0.7,1)+xlim(1,11)+
  labs(x = "Sample size per group", y = "Estimated lambda", title="Delta = 0.15") +
  scale_x_continuous(breaks=seq(0,12,2))+
  geom_hline(yintercept=0.9, linetype="dashed", color = "red")+
  theme_bw()+theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

p3<-ggplot(lambda.data.0.18, aes(x=n0, y=mean)) + 
  geom_errorbar(aes(ymin=low, ymax=up), width=0.2) +
  geom_line(lwd=1, aes(linetype=method)) +
  geom_point()+
  ylim(0.7,1)+xlim(1,11)+
  labs(x = "Sample size per group", y = "Estimated lambda", title="Delta = 0.18") +
  scale_x_continuous(breaks=seq(0,12,2))+
  geom_hline(yintercept=0.9, linetype="dashed", color = "red")+
  theme_bw()+theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

p4<-ggplot(lambda.data.unif, aes(x=n0, y=mean)) + 
  geom_errorbar(aes(ymin=low, ymax=up), width=0.2) +
  geom_line(lwd=1, aes(linetype=method)) +
  geom_point()+
  ylim(0.7,1)+xlim(1,11)+
  labs(x = "Sample size per group", y = "Estimated lambda", title="Delta ~ unif(0.1,0.2)") +
  scale_x_continuous(breaks=seq(0,12,2))+
  geom_hline(yintercept=0.9, linetype="dashed", color = "red")+
  theme_bw()+theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

pdf("results/lambda2by2.pdf")
grid.arrange(p1,p2,p3,p4,ncol=2)
dev.off()