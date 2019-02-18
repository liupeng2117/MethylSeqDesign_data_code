setwd("D:/research/MethySeq/")
source("code/functions.R")
library(limma)
library(qvalue)
library(pi0) #estimate pi0(lambda)
library(ggplot2)
library(doParallel)

load(file="data/CLL/CLL_meth51.rdata") # the raw data is publicly available on GEO, the accession number is GSE66167

############
# differental methylation analysis
############
## Use full data
N.a<-c(20,40,60,80,100,150,200,250)
N.b<-c(4,8,12,16,20,30,40,50)

pilot.N.a<-43
pilot.N.b<-8

target.N.a.ref<-c(pilot.N.a,N.a[N.a>pilot.N.a & N.b>pilot.N.b])
target.N.b.ref<-c(pilot.N.b,N.b[N.a>pilot.N.a & N.b>pilot.N.b])
cases.id<-1:pilot.N.a
controls.id<-44:(44+pilot.N.b-1)
col.coverage<-3*c(cases.id,controls.id)+2
col.numC<-3*c(cases.id,controls.id)+3
coverage<-meth[,col.coverage]
methyl.count<-meth[,col.numC]
coverage<-as.matrix(coverage)
methyl.count<-as.matrix(methyl.count)
  
# Differential methylation analysis
dmr.result<-DMR.analysis(n=c(pilot.N.a, pilot.N.b), cov.matrix = coverage, methyl.matrix = methyl.count)
power.result<-Estimate.EDR.from.pilot(res=dmr.result, thresh.p=0.005,
                                        N0=c(pilot.N.a, pilot.N.b), target.N=cbind(target.N.a.ref, target.N.b.ref),
                                        FDR=0.05, M=20)
edr.ref<-power.result$EDR
fdr.ref<-power.result$FDR


#assign 5 cores
cl<-makeCluster(4)
registerDoParallel(cl)

## predicted power
pilot.N.a.all<-c(10,15,20,30)
pilot.N.b.all<-c(2,3,4,6)
rep<-10
pre.data<-foreach(n = 1:length(pilot.N.a.all), .packages=c("foreach","pi0")) %dopar% {
  
  pilot.N.a<-pilot.N.a.all[n]
  pilot.N.b<-pilot.N.b.all[n]
  
  target.N.a<-c(pilot.N.a,N.a[N.a>pilot.N.a & N.b>pilot.N.b])
  target.N.b<-c(pilot.N.b,N.b[N.a>pilot.N.a & N.b>pilot.N.b])
  
  result.pre<-foreach (times = 1:rep) %do% {
    cases.id<-sample(1:43,pilot.N.a)
    controls.id<-sample(44:51,pilot.N.b)
    col.coverage<-3*c(cases.id,controls.id)+2
    col.numC<-3*c(cases.id,controls.id)+3
    coverage<-meth[,col.coverage]
    methyl.count<-meth[,col.numC]
    coverage<-as.matrix(coverage)
    methyl.count<-as.matrix(methyl.count)
    
    # Differential methylation analysis
    dmr.result<-simulate.methyl.data(n=c(pilot.N.a, pilot.N.b), cov.matrix = coverage, methyl.matrix = methyl.count)
    
    power.result<-Estimate.EDR.from.pilot(res=dmr.result, thresh.p=0.005,
                                          N0=c(pilot.N.a, pilot.N.b), target.N=cbind(target.N.a, target.N.b),
                                          FDR=0.05, M=20)
    y<-rbind(EDR=power.result$EDR, FDR=power.result$FDR, time=b[3])
    return(y)
  }
}
save(pre.data,file="data/CLL/pre51_converge.Rdata")

data.cal<-function(n,result.pre){
  pilot.N.a<-pilot.N.a.all[n]
  pilot.N.b<-pilot.N.b.all[n]
  
  target.N.a<-c(pilot.N.a,N.a[N.a>pilot.N.a & N.b>pilot.N.b])
  target.N.b<-c(pilot.N.b,N.b[N.a>pilot.N.a & N.b>pilot.N.b])
  
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
  edr<-as.factor(c(rep("reference",length(target.N.a.ref)),rep("predicted",length(target.N.a))))
  mean<-c(edr.ref,edr.Zw.inflatN.mean)
  fdrcontrol<-c(fdr.ref,fdr.pre.mean)
  ci.up<-c(edr.ref,edr.Zw.inflatN.mean+1.96*se2)
  ci.dw<-c(edr.ref,edr.Zw.inflatN.mean-1.96*se2)
  pilot.N<-paste(pilot.N.a, "vs", pilot.N.b)
  targetN<-c(paste(target.N.a.ref, "&", target.N.b.ref) , paste(target.N.a, "&", target.N.b))
  controlN<-c(target.N.b.ref,target.N.b)
  sum.data<-data.frame(pilot.N,edr,targetN,controlN,mean,ci.up,ci.dw,fdrcontrol)
  return(sum.data)
}

result<-lapply(1:length(pilot.N.a.all),function(x) data.cal(x,pre.data[[x]]))

levels<-paste(pilot.N.a.all,"vs",pilot.N.b.all)
result.data<-Reduce("rbind",result)
result.data$pilot.N<-factor(result.data$pilot.N,levels=levels)
save(result.data,file="data/CLL/CLL_true_vs_predicted_power51_converge.Rdata")

#plot
pdf("results/CLL_true_vs_predicted_power51_converge.pdf", width = 10, height = 4)
ggplot(result.data, aes(x=controlN, y=mean, colour=edr,shape=edr,linetype=edr)) + 
  geom_errorbar(aes(ymin=ci.up, ymax=ci.dw), width=1) +
  geom_line(lwd=1) +
  geom_point()+
  ylim(0,1)+xlim(0,50)+
  labs(x = "Target N", y = "EDR") +
  scale_linetype_manual(values=c("dotted","solid"))+
  theme_bw()+theme(legend.position = "bottom")+
  facet_grid(~pilot.N,labeller=label_both)
dev.off()

