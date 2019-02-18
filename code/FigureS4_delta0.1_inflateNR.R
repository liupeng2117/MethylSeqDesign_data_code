library(pi0) #estimate pi0(lambda)
library(ggplot2) 
library(VGAM) # estimate fai
library(fitdistrplus) # fit empirical distribution
library(doParallel) # foreach
library(knitr)
library(plyr) #map values of pilotN based on key
library(dplyr)
library(DropletUtils)

setwd("D:/research/MethySeq/")
source("code/functions.R")
load(file="data/Mouse/meth.regional.count.R")
load(file="data/Mouse/empirical.fit.rdata")

pilot.N.all<-c(2,4,6,8,9,10)
N<-c(2,6,10,15,25,50)
rep<-5
G<-10000
n.DE<-1000
#delta<-runif(n.DE,0.1,0.2)
delta<-rep(0.1, n.DE)
prop.all<-c(0.05,0.1,0.2,0.4,0.6,0.8,1.0)
set.seed(123)
dmr<-sample(1:G,n.DE)


#Assign 8 cores
cl<-makeCluster(8)
registerDoParallel(cl)
  
result<-vector("list",length(pilot.N.all)*length(prop.all))

  
for(i in 1:length(prop.all)){
  prop<-prop.all[i]
  
  #True EDR
  result.true<-foreach(j = 1:length(N), .packages = "DropletUtils") %dopar% {
  n<-N[j]
  edr.Zw<-rep(0,rep)
  fdr.Zw<-rep(0,rep)
  dmr.N<-rep(0,rep)
  # simulate dataset with D samples and G CpG regions
  for(times in 1:rep) {
    #simulate data
    simulated.result<-simulate.methyl.data(n=n, G=G, delta=delta, dmr=dmr, prop=prop)
    #DE analysis
    a=proc.time()
    dmr.result<-DMR.analysis(N0=n, cov.matrix=simulated.result$coverage, methyl.matrix=simulated.result$methyl.count, R=1/4, pilot.depth=250/8)
    # Calculate dmr while controling fdr at 0.05
    dmr.Zw<-get.dm.regions(p=dmr.result$p.values, level=0.05, dmr=dmr)
    # Calculate the observed discovery rate
    if(is.na(dmr.Zw[1])){ 
        edr.Zw[times]<-0
        } else {
        edr.Zw[times]<-length(intersect(dmr.Zw,dmr))/length(dmr)
        fdr.Zw[times]<-1-length(intersect(dmr.Zw,dmr))/length(dmr.Zw)
        dmr.N[times]<-length(dmr.Zw)
        result.Zw<-round(rbind(edr.Zw,fdr.Zw,dmr.N),digits=3)
        rownames(result.Zw)<-c("edr","fdr","dmr N")
      }
    }
    return(result.Zw)
  }
  edr.Zw.mean<-unlist(lapply(result.true,function(x) mean(x[1,])))
  edr.Zw.sd<-unlist(lapply(result.true,function(x) sd(x[1,])))
  fdr.Zw.mean<-unlist(lapply(result.true,function(x) mean(x[2,])))
   
  #Estimate EDR
  for(k in 1:length(pilot.N.all)){
      pilot.N<-pilot.N.all[k]
      # power prediction
      result.pre<-foreach (times = 1:rep,.packages = c("pi0", "DropletUtils")) %dopar% {
      #simulate data
      simulated.result<-simulate.methyl.data(n=pilot.N, G=G, delta=delta, dmr=dmr, prop=1)
        
      #DE analysis
      a=proc.time()
      dmr.result<-DMR.analysis(N0=pilot.N,cov.matrix=simulated.result$coverage, methyl.matrix=simulated.result$methyl.count, R=(1/4)*prop, pilot.depth=250/8)
      #estimate EDR
      power.result<-Estimate.EDR.from.pilot(res=dmr.result, thresh.p=0.005,
                                              N0=pilot.N, target.N=N, FDR=0.05, M=10)
      b=proc.time()-a
      y<-rbind(EDR=power.result$EDR, FDR=power.result$FDR, time=b[3])
      return(y)
    }
    
    #Summarize estimation result
    # edr
    result.pre<-result.pre[!unlist(lapply(result.pre,is.null))]
    result.mean<-Reduce("+",result.pre)/length(result.pre)
    result.pre.edr<-t(sapply(result.pre, function(x) x[1,]))
    edr.pre.mean<-result.mean[1,]
    edr.pre.sd<-apply(result.pre.edr,2,sd)
      
    # fdr
    fdr.pre.mean<-result.mean[2,]
    
    #run time
    time.onerun<-mean(result.mean[3,])
      
    # combine to a dataset
    se1<-edr.Zw.sd/sqrt(rep)
    se2<-edr.pre.sd/sqrt(rep)
    edr<-as.factor(c(rep("true",length(N)),rep("predicted",length(N))))
    targetN<-c(N,N)
    mean<-c(edr.Zw.mean,edr.pre.mean)
    fdrcontrol<-c(fdr.Zw.mean,fdr.pre.mean)
    se<-c(se1,se2)
    sum.data<-data.frame(prop=prop,
                         pilot.N=pilot.N, edr=edr,
                         targetN=targetN, mean=mean, 
                         se=se, fdrcontrol=fdrcontrol,
                         time=time.onerun)
    result[[(i-1)*length(pilot.N.all)+k]]<-sum.data
  }
}

new<-c("2 vs 2","4 vs 4","6 vs 6","8 vs 8", "9 vs 9", "10 vs 10")
result.data<-Reduce("rbind",result)
key.values<-data.frame(old=pilot.N.all,new=new)
# change new from factor to character
key.values$new<-as.character(key.values$new)
result.data[, 2] <- mapvalues(result.data[, 2], from = key.values$old, to = key.values$new)
result.data$pilot.N<-factor(result.data$pilot.N,levels=new)


#RMSE
result.true<-result.data %>%
  filter(edr=="true")
result.predicted<-result.data %>%
  filter(edr=="predicted")

result.predicted[,"true"]<-result.true[,"mean"]
result.predicted<-result.predicted %>%
  mutate(squareerror=(mean-true)^2)
result.predicted<-result.predicted %>%
  group_by(prop,pilot.N)%>%
  mutate(rmse=sqrt(sum(squareerror)/length(squareerror)))

result.predicted<-result.predicted %>% 
  dplyr::select(prop,pilot.N,rmse) # select function is covered by MASS
zw<-as.data.frame(result.predicted)
zw<-unique(zw)

save(result, result.data, zw,file="data/simulation/powerplotrawdataCBUM+CDD0.005_20iter_cleaned_delta0.1_infateNR.Rdata")

#plot
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols<-gg_color_hue(2)
pdf("results/CBUM+CDD0.005_20iter_delta0.1_inflateNR.pdf", width = 9, height = 9)
ggplot(result.data, aes(x=targetN, y=mean, colour=edr,shape=edr,linetype=edr)) + 
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), width=1) +
  geom_line(lwd=1) +
  geom_point()+
  ylim(0,1)+xlim(0,50)+
  labs(x = "Target N", y = "EDR") +
  scale_linetype_manual(values=c("dotted","solid"))+
  theme_bw()+theme(legend.position = "bottom")+
  facet_grid(prop~pilot.N,labeller=label_both)
dev.off()


