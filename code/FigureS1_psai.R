library(pi0) #estimate pi0(lambda)
library(ggplot2) 
library(VGAM) # estimate fai
library(fitdistrplus) # fit empirical distribution
library(doParallel) # foreach
library(knitr)
library(plyr) #map values of pilotN based on key
library(reshape2) # reshape
library(gridExtra)
library(DropletUtils) #subsampling

setwd("D:/research/MethySeq/")
load(file="data/Mouse/meth.regional.count.R")

##=================================##
#=====  Sample size=12      ========#
##=================================##
set.seed(123)
G=1
D=12
fai=0.04838581
seq.coverage=seq(5,32,by=3)
seq.numC=seq(6,33,by=3)
coverage<-meth[,seq.coverage]
size.mean<-6.286579

# First calculate sequencing depth and theta-hat
#Then downsample our sequencing depth to 5%, 10%, 20%,..., 80%
folds<-c(0.05,0.1,0.2,0.4,0.6,0.8,1)
ratio.new<-matrix(0,G,length(folds))
ratio.new.dn<-matrix(0,G,length(folds))
ratio.new.up<-matrix(0,G,length(folds))
raw.ratio<-list()
for(i in 1:length(folds)){
  fold=folds[i]
  # calculate new seq.depth and new coverage level of each region
  ratio40<-matrix(0,G,40)
  for(k in 1:40){
    sim.coverage<-matrix(0,G,D)
    for(h in 1:G){
      sim.coverage[h,]<-rnbinom(D,mu=median(as.matrix(coverage)),size=size.mean)
    }
    
    median(sim.coverage)
    
    # generate coverage from multinomial distribution
    coverage.new<-downsampleMatrix(sim.coverage, prop = fold)
    #coverage.new<-matrix(rmultinom(n=1,size=seq.depth.new,prob=theta.hat.vector),ncol=D)
    #calculate the updated value
    data.new<-matrix(0,G,D)
    for(j in 1:D){
      m<-coverage.new[,j]
      data.new[,j]<-m/(1+(m-1)*fai)
    }
    
    mean.group1<-apply(data.new[,1:(D/2),drop=F],1,mean)
    mean.group2<-apply(data.new[,(D/2+1):D,drop=F],1,mean)
    ratio40[,k]<-(mean.group1+mean.group2)/(mean.group1*mean.group2)
  }
  ratio.median<-apply(ratio40,1,median)
  ratio.q<-apply(ratio40,1,function(x) quantile(x,probs=c(0.25,0.75)))
  raw.ratio[[i]]<-ratio40
  ratio.new[,i]<-ratio.median
  ratio.new.dn[,i]<-ratio.q[1,]
  ratio.new.up[,i]<-ratio.q[2,]
}

save(ratio.new,ratio.new.up,ratio.new.dn,raw.ratio,file="samples12.Rdata")

cov.g<-apply(sim.coverage,1,mean)
id<-which.min(abs(cov.g-median(cov.g)))
a<-sapply(raw.ratio,function(x) x[id,])

df<-cbind(folds=folds, ratio=as.data.frame(t(a)))
#df.long<-melt(df,id.var="folds")
mean<-apply(df[2:ncol(df)],1,mean)
se<-apply(df[2:ncol(df)],1,sd)
ci.up<-mean+1.96*se
ci.dn<-mean-1.96*se
df.wide<-data.frame(folds=df$folds,mean,se,ci.up,ci.dn)

##=================================##
#=====  Sample size=6     ========#
##=================================##
set.seed(123)
D=6

folds<-c(0.05,0.1,0.2,0.4,0.6,0.8,1)
ratio.new<-matrix(0,G,length(folds))
ratio.new.dn<-matrix(0,G,length(folds))
ratio.new.up<-matrix(0,G,length(folds))
raw.ratio<-list()
for(i in 1:length(folds)){
  fold=folds[i]
  # calculate new seq.depth and new coverage level of each region
  ratio40<-matrix(0,G,40)
  for(k in 1:40){
    sim.coverage<-matrix(0,G,D)
    for(h in 1:G){
      sim.coverage[h,]<-rnbinom(D,mu=median(as.matrix(coverage)),size=size.mean)
    }
    
    median(sim.coverage)
    
    # generate coverage from multinomial distribution
    coverage.new<-downsampleMatrix(sim.coverage, prop = fold)
    #coverage.new<-matrix(rmultinom(n=1,size=seq.depth.new,prob=theta.hat.vector),ncol=D)
    #calculate the updated value
    data.new<-matrix(0,G,D)
    for(j in 1:D){
      m<-coverage.new[,j]
      data.new[,j]<-m/(1+(m-1)*fai)
    }
    
    mean.group1<-apply(data.new[,1:(D/2),drop=F],1,mean)
    mean.group2<-apply(data.new[,(D/2+1):D,drop=F],1,mean)
    ratio40[,k]<-(mean.group1+mean.group2)/(mean.group1*mean.group2)
  }
  ratio.median<-apply(ratio40,1,median)
  ratio.q<-apply(ratio40,1,function(x) quantile(x,probs=c(0.25,0.75)))
  raw.ratio[[i]]<-ratio40
  ratio.new[,i]<-ratio.median
  ratio.new.dn[,i]<-ratio.q[1,]
  ratio.new.up[,i]<-ratio.q[2,]
}

save(ratio.new,ratio.new.up,ratio.new.dn,raw.ratio,file="samples3.Rdata")

cov.g<-apply(sim.coverage,1,mean)
id<-which.min(abs(cov.g-median(cov.g)))
a<-sapply(raw.ratio,function(x) x[id,])

df0<-cbind(folds=folds, ratio=as.data.frame(t(a)))
#df0.long<-melt(df0,id.var="folds")
mean<-apply(df0[2:ncol(df0)],1,mean)
se<-apply(df0[2:ncol(df0)],1,sd)
ci.up<-mean+1.96*se
ci.dn<-mean-1.96*se
df0.wide<-data.frame(folds=df0$folds,mean,se,ci.up,ci.dn)

##=================================##
#=====  Sample size=18      ========#
##=================================##
set.seed(123)
D=18

folds<-c(0.05,0.1,0.2,0.4,0.6,0.8,1)
ratio.new<-matrix(0,G,length(folds))
ratio.new.dn<-matrix(0,G,length(folds))
ratio.new.up<-matrix(0,G,length(folds))
raw.ratio<-list()
for(i in 1:length(folds)){
  fold=folds[i]
  # calculate new seq.depth and new coverage level of each region
  ratio40<-matrix(0,G,40)
  for(k in 1:40){
    sim.coverage<-matrix(0,G,D)
    for(h in 1:G){
      sim.coverage[h,]<-rnbinom(D,mu=median(as.matrix(coverage)),size=size.mean)
    }
    
    median(sim.coverage)
    
    # generate coverage from multinomial distribution
    coverage.new<-downsampleMatrix(sim.coverage, prop = fold)
    #coverage.new<-matrix(rmultinom(n=1,size=seq.depth.new,prob=theta.hat.vector),ncol=D)
    #calculate the updated value
    data.new<-matrix(0,G,D)
    for(j in 1:D){
      m<-coverage.new[,j]
      data.new[,j]<-m/(1+(m-1)*fai)
    }
    
    mean.group1<-apply(data.new[,1:(D/2),drop=F],1,mean)
    mean.group2<-apply(data.new[,(D/2+1):D,drop=F],1,mean)
    ratio40[,k]<-(mean.group1+mean.group2)/(mean.group1*mean.group2)
  }
  ratio.median<-apply(ratio40,1,median)
  ratio.q<-apply(ratio40,1,function(x) quantile(x,probs=c(0.25,0.75)))
  raw.ratio[[i]]<-ratio40
  ratio.new[,i]<-ratio.median
  ratio.new.dn[,i]<-ratio.q[1,]
  ratio.new.up[,i]<-ratio.q[2,]
}

save(ratio.new,ratio.new.up,ratio.new.dn,raw.ratio,file="samples9.Rdata")

cov.g<-apply(sim.coverage,1,mean)
id<-which.min(abs(cov.g-median(cov.g)))
a<-sapply(raw.ratio,function(x) x[id,])

df2<-cbind(folds=folds, ratio=as.data.frame(t(a)))
#df2.long<-melt(df2,id.var="folds")
mean<-apply(df2[2:ncol(df2)],1,mean)
se<-apply(df2[2:ncol(df2)],1,sd)
ci.up<-mean+1.96*se
ci.dn<-mean-1.96*se
df2.wide<-data.frame(folds=df2$folds,mean,se,ci.up,ci.dn)

##=================================##
#=====  Sample size=24     ========#
##=================================##
set.seed(123)
D=24
folds<-c(0.05,0.1,0.2,0.4,0.6,0.8,1)
ratio.new<-matrix(0,G,length(folds))
ratio.new.dn<-matrix(0,G,length(folds))
ratio.new.up<-matrix(0,G,length(folds))
raw.ratio<-list()
for(i in 1:length(folds)){
  fold=folds[i]
  # calculate new seq.depth and new coverage level of each region
  ratio40<-matrix(0,G,40)
  for(k in 1:40){
    sim.coverage<-matrix(0,G,D)
    for(h in 1:G){
      sim.coverage[h,]<-rnbinom(D,mu=median(as.matrix(coverage)),size=size.mean)
    }
    
    median(sim.coverage)
    
    # generate coverage from multinomial distribution
    coverage.new<-downsampleMatrix(sim.coverage, prop = fold)
    #coverage.new<-matrix(rmultinom(n=1,size=seq.depth.new,prob=theta.hat.vector),ncol=D)
    #calculate the updated value
    data.new<-matrix(0,G,D)
    for(j in 1:D){
      m<-coverage.new[,j]
      data.new[,j]<-m/(1+(m-1)*fai)
    }
    
    mean.group1<-apply(data.new[,1:(D/2),drop=F],1,mean)
    mean.group2<-apply(data.new[,(D/2+1):D,drop=F],1,mean)
    ratio40[,k]<-(mean.group1+mean.group2)/(mean.group1*mean.group2)
  }
  ratio.median<-apply(ratio40,1,median)
  ratio.q<-apply(ratio40,1,function(x) quantile(x,probs=c(0.25,0.75)))
  raw.ratio[[i]]<-ratio40
  ratio.new[,i]<-ratio.median
  ratio.new.dn[,i]<-ratio.q[1,]
  ratio.new.up[,i]<-ratio.q[2,]
}

save(ratio.new,ratio.new.up,ratio.new.dn,raw.ratio,file="samples12.Rdata")

cov.g<-apply(sim.coverage,1,mean)
id<-which.min(abs(cov.g-median(cov.g)))
a<-sapply(raw.ratio,function(x) x[id,])

df3<-cbind(folds=folds, ratio=as.data.frame(t(a)))
#df3.long<-melt(df3,id.var="folds")
mean<-apply(df3[2:ncol(df3)],1,mean)
se<-apply(df3[2:ncol(df3)],1,sd)
ci.up<-mean+1.96*se
ci.dn<-mean-1.96*se
df3.wide<-data.frame(folds=df3$folds,mean,se,ci.up,ci.dn)

save(df0,df0.wide,df,df.wide,df2,df2.wide,df3,df3.wide,file="data/simulation/Stability_quantity.Rdata")

###======================###
# plot #
###======================###

plot1<-ggplot(df0.wide, aes(x=folds, y=mean)) + 
  geom_errorbar(aes(ymin=ci.up, ymax=ci.dn), width=.05) +
  geom_line() +
  geom_point() +
  theme_bw() + ylim(0,0.3)+ theme(plot.title = element_text(hjust = 0.5))+
  xlab(expression("R"[j]*"/R"[0])) + ylab(expression(Psi)) + ggtitle("3 vs. 3")


plot2<-ggplot(df.wide, aes(x=folds, y=mean)) + 
  geom_errorbar(aes(ymin=ci.up, ymax=ci.dn), width=.05) +
  geom_line() +
  geom_point() +
  theme_bw() + ylim(0,0.3)+ theme(plot.title = element_text(hjust = 0.5))+
  xlab(expression("R"[j]*"/R"[0])) + ylab(expression(Psi)) + ggtitle("6 vs. 6")

plot3<-ggplot(df2.wide, aes(x=folds, y=mean)) + 
  geom_errorbar(aes(ymin=ci.up, ymax=ci.dn), width=.05) +
  geom_line() +
  geom_point() +
  theme_bw() + ylim(0,0.3)+ theme(plot.title = element_text(hjust = 0.5))+
  xlab(expression("R"[j]*"/R"[0])) + ylab(expression(Psi)) + ggtitle("9 vs. 9")

plot4<-ggplot(df3.wide, aes(x=folds, y=mean)) + 
  geom_errorbar(aes(ymin=ci.up, ymax=ci.dn), width=.05) +
  geom_line() +
  geom_point() +
  theme_bw() + ylim(0,0.3)+ theme(plot.title = element_text(hjust = 0.5))+
  xlab(expression("R"[j]*"/R"[0])) + ylab(expression(Psi)) + ggtitle("12 vs. 12")

pdf("results/stability_quantity2.pdf",width = 6,height = 4)
grid.arrange(plot1, plot2, plot3, plot4, ncol=2)
dev.off()
