setwd("D:/research/MethySeq/")
source("code/functions.R")
library(limma)
library(qvalue)
library(pi0) #estimate pi0(lambda)
library(ggplot2)
library(doParallel)
load(file="data/Mouse/meth.regional.count.R")

############
# differental methylation analysis
############
## Use full data
N<-c(4,12,20,30,50,100)

#total 10 samples. V5 and P12 were removed from furthur analysis due to data quality issues.
pilot.N<-5

target.N.ref<-c(pilot.N,N[N>pilot.N])
seq.coverage<-seq(5,32,by=3)
seq.numC<-seq(6,33,by=3)
coverage<-meth[,seq.coverage]
methyl.count<-meth[,seq.numC]
coverage<-as.matrix(coverage)
methyl.count<-as.matrix(methyl.count)
  
# Differential methylation analysis
dmr.result<-DMR.analysis(n=pilot.N, cov.matrix = coverage, methyl.matrix = methyl.count)
power.result<-Estimate.EDR.from.pilot(res=dmr.result, thresh.p=0.005,
                                        N0=pilot.N, target.N=target.N.ref,
                                        FDR=0.05, M=20)

edr.ref<-power.result$EDR
fdr.ref<-power.result$FDR

#Use subsampled pilot data
#assign 5 cores
cl<-makeCluster(5)
registerDoParallel(cl)

## predicted power
pilot.N.all<-c(2,3,4,5)

rep<-10
pre.data<-foreach(n = 1:length(pilot.N.all), .packages=c("foreach","pi0")) %dopar% {
  
  pilot.N<-pilot.N.all[n]
  target.N<-c(pilot.N,N[N>pilot.N])
  
  result.pre<-foreach (times = 1:rep) %do% {
    cases.id<-sample(1:5,pilot.N)
    controls.id<-sample(6:10,pilot.N)
    col.coverage<-3*c(cases.id,controls.id)+2
    col.numC<-3*c(cases.id,controls.id)+3
    coverage<-meth[,col.coverage]
    methyl.count<-meth[,col.numC]
    coverage<-as.matrix(coverage)
    methyl.count<-as.matrix(methyl.count)
    
    # Differential methylation analysis
    dmr.result<-DMR.analysis(n=pilot.N, cov.matrix = coverage, methyl.matrix = methyl.count)
    power.result<-Estimate.EDR.from.pilot(res=dmr.result, thresh.p=0.005,
                                          N0=pilot.N, target.N=target.N,
                                          FDR=0.05, M=20)
    y<-rbind(EDR=power.result$EDR, FDR=power.result$FDR)
    return(y)
  }
}
save(pre.data,file="data/Mouse/Mouse_pre.Rdata")

data.cal<-function(n,result.pre){
  pilot.N<-pilot.N.all[n]
  target.N<-c(pilot.N,N[N>pilot.N & N>pilot.N])
  
  # edr
  result.pre<-result.pre[!unlist(lapply(result.pre,is.null))]
  result.mean<-Reduce("+",result.pre)/length(result.pre)
  result.pre.edr<-t(sapply(result.pre, function(x) x[1,]))
  edr.Zw.inflatN.mean<-result.mean[1,]
  edr.Zw.inflatN.sd<-apply(result.pre.edr,2,sd)
  
  # fdr
  fdr.pre.mean<-result.mean[2,]
  
  #dataset
  #se1<-edr.Zw.sd/sqrt(rep)
  se2<-edr.Zw.inflatN.sd/sqrt(rep)
  edr<-as.factor(c(rep("reference",length(target.N.ref)),rep("predicted",length(target.N))))
  mean<-c(edr.ref,edr.Zw.inflatN.mean)
  fdrcontrol<-c(fdr.ref,fdr.pre.mean)
  ci.up<-c(edr.ref,edr.Zw.inflatN.mean+1.96*se2)
  ci.dw<-c(edr.ref,edr.Zw.inflatN.mean-1.96*se2)
  pilot.N<-paste(pilot.N, "vs", pilot.N)
  targetN<-c(paste(target.N.ref, "&", target.N.ref) , paste(target.N, "&", target.N))
  controlN<-c(target.N.ref,target.N)
  sum.data<-data.frame(pilot.N,edr,targetN,controlN,mean,ci.up,ci.dw,fdrcontrol)
  return(sum.data)
}

result<-lapply(1:length(pilot.N.all),function(x) data.cal(x,pre.data[[x]]))

levels<-paste(pilot.N.all,"vs",pilot.N.all)
result.data<-Reduce("rbind",result)
result.data$pilot.N<-factor(result.data$pilot.N,levels=levels)
save(result.data,file="data/Mouse/Mouse_ref_vs_pred.Rdata")

#plot
pdf("results/Mouse_ref_vs_pre.pdf", width = 8, height = 4)
ggplot(result.data, aes(x=controlN, y=mean, colour=edr,shape=edr,linetype=edr)) + 
  geom_errorbar(aes(ymin=ci.up, ymax=ci.dw), width=1) +
  geom_line(lwd=1) +
  geom_point()+
  ylim(0,1)+xlim(0,50)+
  labs(x = "Target N", y = "EDR") +
  scale_linetype_manual(values=c("dotted","solid"))+
  scale_shape_manual(name = "EDR",labels = c("Predicted", "Reference"),values=c(16,17))+  
  theme_bw()+theme(legend.position = "bottom")+
  facet_grid(~pilot.N,labeller=label_both)
dev.off()

