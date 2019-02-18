##############
# Functions when generating the datasets
##############
##### simulate data #######

simulate.methyl.data<-function(n,G,delta,low.bound=0,up.bound=1,dmr,prop=1){ # delta follows uniform distribution
  D=2*n
  ### Simulate coverage
  #draw G mean coverage miu from empirical distribution
  get.c<-function(G){emp.c <- as.numeric(quantile(fitemp.c, runif(G)))}
  c.sim<-get.c(G)
  c.sim=prop*c.sim
  
  # simulate coverage using mean drawing from empirical distribution and dispersion estimate from data
  #sim.coverage<-matrix(sample(coverage,G*D,replace = T),G,D)
  sim.coverage<-matrix(0,G,D)
  for(i in 1:G){
    sim.coverage[i,]<-rnbinom(D,mu=c.sim[i],size=size.mean)
  }
  
  # draw G mean p
  get.p<-function(G){emp.p <- as.numeric(quantile(fitemp.p, runif(G,fitemp.p(low.bound),fitemp.p(up.bound))))}
  p.sim<-get.p(G)
  
  # If dmr, then methylation difference is delta
  # change the pvalue of cases
  p.sim.cases<-rep(0,G)
  for(i in 1:G){
    if(i %in% dmr){
      pos<-which(i==dmr)
      p.sim.cases[i]<-ifelse(p.sim[i]>=0.5,p.sim[i]-delta[pos],p.sim[i]+delta[pos])
    } else {
      p.sim.cases[i]<-p.sim[i]
    }
  }
  
  #calculate alpha and beta
  theta<-1/fai.mean-1
  alpha.cases=p.sim.cases*theta
  beta.cases=theta-alpha.cases
  alpha.controls<-p.sim*theta
  beta.controls<-theta-alpha.controls
  
  # simulate p for each cell
  p.sim.all<-matrix(0,G,D)
  for (j in 1:G) {
    if(j %in% dmr){
      p.sim.all[j,1:(D/2)]<-rbeta(D/2,alpha.cases[j],beta.cases[j])
      p.sim.all[j,(D/2+1):D]<-rbeta(D/2,alpha.controls[j],beta.controls[j])
    } else {
      p.sim.all[j,]<-rbeta(D,alpha.controls[j],beta.controls[j])
    }
  }
  
  ### generagte methylation level x based on p and coverage
  x<-matrix(0,G,D)
  for(i in 1:D){
    for(j in 1:G){
      x[j,i]<-rbinom(1,sim.coverage[j,i],p.sim.all[j,i])
    }
  }
  return(list(coverage=sim.coverage,methyl.count=x,p.controls=p.sim,p.cases=p.sim.cases))
}

##############################################
# Comparing different tests
##############################################
# Beat + t-test
B.ttest<-function(sim.coverage,sim.meth,a=1,D=10){
  beta<-(sim.meth+a)/(sim.coverage+2*a)
  G<-nrow(sim.coverage)
  p.Bt<-rep(0,G)
  for(i in 1:G){
    if(length(unique(beta[i,]))<=2){
      a<-runif(D)/10^(4)
      p.Bt[i]<-t.test(beta[i,1:(D/2)]+a[1:(D/2)],beta[i,(D/2+1):D]+a[(D/2+1):D],var.equal = T)$p.value
    } else{
    p.Bt[i]<-t.test(beta[i,1:D/2],beta[i,(D/2+1):D],var.equal = T)$p.value
    }
  }
  return(p.Bt)
}
# M-value + t-test
M.ttest<-function(sim.coverage,sim.meth,a=1,D=10){
  M.value<-log2((sim.meth+a)/(sim.coverage-sim.meth+a))
  G<-nrow(sim.coverage)
  p.Mt<-rep(0,G)
  for(i in 1:G){
    if(length(unique(M.value[i,]))<=2){ # in case that all values of t test are the same
      a<-runif(D)/10^(4)
      p.Mt[i]<-t.test(M.value[i,1:(D/2)]+a[1:(D/2)],M.value[i,(D/2+1):D]+a[(D/2+1):D],var.equal = T)$p.value
    } else{
    p.Mt[i]<-t.test(M.value[i,1:D/2],M.value[i,(D/2+1):D],var.equal = T)$p.value
    }
  }
  return(p.Mt)
}


# Z-value + t-test
Z.ttest<-function(sim.coverage,sim.meth,a=1,D=10){
  Z.value<-asin(2*(sim.meth+a)/(sim.coverage+2*a)-1)
  G<-nrow(sim.coverage)
  p.Zt<-rep(0,G)
  for(i in 1:G){
    if(length(unique(Z.value[i,]))<=2){ # in case that all values of t test are the same
      a<-runif(D)/10^(4)
      p.Zt[i]<-t.test(Z.value[i,1:(D/2)]+a[1:(D/2)],Z.value[i,(D/2+1):D]+a[(D/2+1):D],var.equal = T)$p.value
    } else{
      p.Zt[i]<-t.test(Z.value[i,1:(D/2)],Z.value[i,(D/2+1):D],var.equal = T)$p.value
    }
  }
  return(p.Zt)
}


# Z-value + wald-test
# Calculate the value

# Z-value + wald-test
Z.wald <- function(x, cov.matrix, sim.meth, fai, a = 1) {
  D <- sum(x)
  G <- dim(cov.matrix)[1]
  if (dim(cov.matrix)[2] != D) {
    stop("sample size and coverage data matrix columns differ")
  }
  if (dim(sim.meth)[2] != D) {
    stop("sample size and methylation data matrix columns differ")
  }
  n.groups <- length(x)
  if (n.groups > 2) {
    stop("Currently more than 2 groups is not allowed")
  }
  
  # Assign each group size
  for (i in 1:n.groups) {
    assign(paste0("D", i), x[i])
  }
  
  ### calculate beta1-hat
  chi2 <- rep(0, G)
  beta0.hat<-rep(0,G)
  beta1.hat<-rep(0,G)
  
  for (i in 1:G) {
    tryCatch({
      X <- matrix(c(rep(1, D), rep(1, D1), rep(0, D2)), ncol = 2)
      V.inverse <- diag(cov.matrix[i, ]/(1 + (cov.matrix[i, ] - 1) * fai))
      Z <- asin(2 * (sim.meth[i, ] + a)/(cov.matrix[i, ] + 2 * a) - 1)
      sigma <- solve(t(X) %*% V.inverse %*% X)
      beta.hat <- sigma %*% t(X) %*% V.inverse %*% Z
      beta0.hat[i]<-beta.hat[1]
      beta1.hat[i]<-beta.hat[2]
      
      C <- diag(rep(1, n.groups))[, -1]
      var <- t(C) %*% sigma %*% C
      chi2[i] <- t(t(C) %*% beta.hat) %*% solve(var) %*% (t(C) %*% beta.hat)
    }, error = function(err) {
      # error handler picks up where error was generated
      print(paste("ERROR:  ", err))
    })
  }
  # Calculate test statistic, p and q values
  p.Zw <- pchisq(chi2, df = n.groups - 1, lower.tail = F)
  return(list(p.Zw = p.Zw, beta0=beta0.hat, beta1=beta1.hat, z.statistic = chi2))
}


## Estimate fai using park's method
estimate.fai <- function(cov.matrix, methyl.matrix, a = 1, x) {
  G <- dim(cov.matrix)[1]
  D <- dim(cov.matrix)[2]
  n.groups <- length(x)
  # Assign each group size
  for (i in 1:n.groups) {
    assign(paste0("D", i), x[i])
  }
  
  if (D == 2) 
    fai.est = 0.001 else {
      fai <- rep(0, G)
      for (i in 1:G) {
        tryCatch({
          V0.inverse <- diag(cov.matrix[i, ])
          # X matrix
          X <- matrix(c(rep(1, D), rep(1, D1), rep(0, D - D1)), ncol = 2)
          if (n.groups > 2) {
            c <- D1
            for (j in 2:(n.groups - 1)) {
              col <- c(rep(0, c), rep(1, get(paste0("D", j))), rep(0, D - c - get(paste0("D", j))))
              X <- cbind(X, col)
              c = c + get(paste0("D", j))
            }
          }
          Z <- asin(2 * (methyl.matrix[i, ] + a)/(cov.matrix[i, ] + 2 * a) - 1)
          sigma <- solve(t(X) %*% V0.inverse %*% X)
          beta.hat0 <- sigma %*% t(X) %*% V0.inverse %*% as.matrix(Z)
          chi2 <- sum(cov.matrix[i, ] * (Z - X %*% beta.hat0)^2)
          sigma2.hat <- chi2/(D - n.groups)
          fai[i] <- (D * (sigma2.hat - 1))/sum(cov.matrix[i, ] - 1)
        }, error = function(err) {
          # error handler picks up where error was generated
          print(paste("MY_ERROR:  ", err))
        })
      }
      fai <- fai[complete.cases(fai)]
      fai <- ifelse(fai < 0, 10^(-6), fai)
      fai <- ifelse(fai > 1, 1 - 10^(-6), fai)
      fai.est <- mean(fai)
    }
  return(fai.est)
}

#=======================
# DMR analysis
#=======================
DMR.analysis <- function(N0, cov.matrix, methyl.matrix, R, pilot.depth, depth_per_lane=250, align_rate=0.5) {
  if (class(N0) != "numeric" | length(N0) > 2) {
    stop("Argument N0 is not correctly specified")
  } else if (length(N0) == 1) {
    N0 = rep(N0, 2)
  }
  if (class(cov.matrix) == "data.frame") {
    cov.matrix = as.matrix(cov.matrix)
  }
  if (class(methyl.matrix) == "data.frame") {
    methyl.matrix = as.matrix(methyl.matrix)
  }
  if (class(cov.matrix) != "matrix") {
    stop("input coverage matrix should be either matrix or data frame")
  }
  if (class(methyl.matrix) != c("matrix")) {
    stop("input coverage matrix should be either matrix or data frame")
  }
  if (!"DropletUtils" %in% rownames(installed.packages())) {
    install.packages("DropletUtils")
  }
  library(DropletUtils)
  print(paste(N0, sep = "", collapse = " vs "))
  pilot.R = (pilot.depth/align_rate)/depth_per_lane
  R = R[R<=pilot.R]
  fai.est <- estimate.fai(x = N0, cov.matrix, methyl.matrix)
  test.result <- Z.wald(x = N0, cov.matrix, methyl.matrix, fai = fai.est, a = 1)
  p.values <- test.result[[1]]
  # calculate ratio of psai
  factorR <- matrix(, nrow = length(p.values), ncol = length(R))
  colnames(factorR) <- R
  for (i in 1:length(R)) {
    cov.matrix1 <- downsampleMatrix(cov.matrix, prop = R[i]/pilot.R)
    data0 = data1 = cov.matrix
    for (j in 1:sum(N0)) {
      m0 <- cov.matrix[, j]
      m1 <- cov.matrix1[, j]
      data0[, j] <- m0/(1 + (m0 - 1) * fai.est)
      data1[, j] <- m1/(1 + (m1 - 1) * fai.est)
    }
    
    mean.group1 <- apply(data0[, 1:N0[1], drop = F], 1, mean)
    mean.group2 <- apply(data0[, (N0[1] + 1):sum(N0), drop = F], 1, mean)
    psai0 <- (mean.group1 + mean.group2)/(mean.group1 * mean.group2)
    
    mean.group1 <- apply(data1[, 1:N0[1], drop = F], 1, mean)
    mean.group2 <- apply(data1[, (N0[1] + 1):sum(N0), drop = F], 1, mean)
    psai1 <- (mean.group1 + mean.group2)/(mean.group1 * mean.group2)
    
    factorR[, i] = psai0/psai1
  }
  
  model <- cbind(beta0 = test.result[[2]], beta1 = test.result[[3]], fai = fai.est)
  return(list(p.values = p.values, model = model, factorR = factorR))
}

#############################
#   Mixture model fitting   #
#############################

# Use censored beta-uniform mixture model to estimate pi0(lambda)
MixtureModel.Fittting.pi0 <- function(p.values, restrict = TRUE, l.upper = 0.95, thresh.p = 0.05) {
  if (!"limma" %in% rownames(installed.packages())) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("limma")
  }
  if (!"pi0" %in% rownames(installed.packages())) {
    devtools::install_github("gitlongor/pi0")
  }
  
  library(limma)
  library(pi0)
  
  fBUM <- function(z, d, lambda) {
    p.values = d
    r = z[1]
    s = z[2]
    Target.Func = -sum(log(lambda + (1 - lambda) * p.values^(r - 1) * (1 - p.values)^(s - 1)/beta(r, s)))
    return(Target.Func)
  }
  
  ## inital value(MME and non-parametric)
  mean.p = mean(p.values, na.rm = T)
  var.p = var(p.values, na.rm = T)
  init.r = ((1 - mean.p) * mean.p^2 - mean.p * var.p)/var.p
  init.s = ((1 - mean.p)^2 * mean.p - (1 - mean.p) * var.p)/var.p
  init = c(max(0, min(init.r, 0.9)), max(1, init.s))
  
  lambda0 = CBUM(p.values, thresh.censor = thresh.p, niter = 1000)
  
  if (attributes(lambda0)$converged == FALSE) {
    print("Fails to estimate lambda, CBUM model does not converge. Use CDD instead")
    lambda0 = convest(p.values)[[1]]
  }
  
  if (restrict) {
    lambda = max(min(lambda0, l.upper), 0.8)
    r.upper = 0.9
    s.upper = Inf
    r.lower = 0
    s.lower = 1
  } else {
    lambda = min(lambda0, 0.99)
    r.upper = 1
    s.upper = Inf
    r.lower = 0
    s.lower = 1
  }
  tryCatch({
    pars = optim(init, fn = fBUM, d = p.values, lambda = lambda, method = "L-BFGS-B", upper = c(r.upper, s.upper), lower = c(r.lower, s.lower))$par
    LL = sum(log(lambda + (1 - lambda) * p.values^(pars[1] - 1) * (1 - p.values)^(pars[2] - 1)/beta(pars[1], pars[2])))
    out.par <- c(lambda, pars, LL)
    names(out.par) <- c("lambda", "r", "s", "LL")
    return(out.par)
  }, error = function(e) {
    print("Fail to optimize r and s, return initial values.")
    LL = sum(log(lambda + (1 - lambda) * p.values^(init[1] - 1) * (1 - p.values)^(init[2] - 1)/beta(init[1], init[2])))
    out.par <- c(lambda, init, LL)
    names(out.par) <- c("lambda", "r", "s", "LL")
    return(out.par)
  })
}

# Use beta-uniform mixture model to estimate pi0(lambda)
MixtureModel.Fittting.pi0.MLE <- function(p.values, restrict = TRUE, l.upper = 0.95, thresh.p = 0.05) {
  if (!"limma" %in% rownames(installed.packages())) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("limma")
  }
  if (!"pi0" %in% rownames(installed.packages())) {
    devtools::install_github("gitlongor/pi0")
  }
  
  library(limma)
  library(pi0)
  
  fBUM <- function(z, d) {
    p.values = d
    lambda = z[1]
    r = z[2]
    s = z[3]
    Target.Func = -sum(log(lambda + (1 - lambda) * p.values^(r - 1) * (1 - p.values)^(s - 1)/beta(r, s)))
    return(Target.Func)
  }
  
  ## inital value(MME and non-parametric)
  mean.p = mean(p.values, na.rm = T)
  var.p = var(p.values, na.rm = T)
  init.r = ((1 - mean.p) * mean.p^2 - mean.p * var.p)/var.p
  init.s = ((1 - mean.p)^2 * mean.p - (1 - mean.p) * var.p)/var.p
  lambda0 = CBUM(p.values, thresh.censor = thresh.p, niter = 1000)
  init = c(min(0.5,max(lambda0,0.99)), max(0, min(init.r, 0.9)), max(1, init.s))
  #if (attributes(lambda0)$converged == FALSE) {
  #  print("Fails to estimate lambda, CBUM model does not converge. Use CDD instead")
  #  lambda0 = convest(p.values)[[1]]
  #}
  
  if (restrict) {
    l.upper = 0.99
    l.lower = 0.8
    r.upper = 0.9
    s.upper = Inf
    r.lower = 0
    s.lower = 1
  } else {
    l.upper = 0.99
    l.lower = 0.5
    r.upper = 1
    s.upper = Inf
    r.lower = 0
    s.lower = 1
  }
  tryCatch({
    pars = optim(init, fn = fBUM, d = p.values, method = "L-BFGS-B", upper = c(l.upper, r.upper, s.upper), lower = c(l.lower, r.lower, s.lower))$par
    LL = sum(log(pars[1] + (1 - pars[1]) * p.values^(pars[2] - 1) * (1 - p.values)^(pars[3] - 1)/beta(pars[2], pars[3])))
    out.par <- c(pars, LL)
    names(out.par) <- c("lambda", "r", "s", "LL")
    return(out.par)
  }, error = function(e) {
    print("Fail to optimize r and s, return initial values.")
    LL = sum(log(init[1] + (1 - init[1]) * p.values^(init[2] - 1) * (1 - p.values)^(init[3] - 1)/beta(init[2], init[3])))
    out.par <- c(init, LL)
    names(out.par) <- c("lambda", "r", "s", "LL")
    return(out.par)
  })
}


# Get DMR when controling true FDR at 0.05 level
get.dm.regions <- function(p, level = 0.05, dmr) {
  p.de <- p[dmr]
  p.ordered <- p[order(p)]
  tdr <- rep(0, length(p))
  for (i in 1:length(p)) {
    tdr[i] = sum(p.de <= p.ordered[i])/sum(p.ordered <= p.ordered[i])
  }
  id <- which(tdr >= 1 - level)
  if (length(id) == 0) {
    return(print("0 DMR found"))
  } else num.de <- max(id)
  claim.dmr <- order(p)[1:num.de]
  return(claim.dmr)
}
##############################
# Estimate lambda from pilot #
##############################
Estimate.lambda.from.pilot <- function(p.values, N0, target.N, thresh.p = 0.005, FDR = 0.05, M = 10) {
  
  if (all(is.na(p.values))) {
    stop("all p-values were NA; nothing to compute")
  }
  if (any(is.na(p.values))) {
    p.values <- p.values[!is.na(p.values)]
  }
  if (min(p.values) == 0) {
    min.nonzero <- min(p.values[p.values > 0])
    p.values[p.values == 0] <- min.nonzero/2
  }
  if (max(p.values) == 1) {
    max.nonone <- max(p.values[p.values < 1])
    p.values[p.values == 1] <- 1 - (1 - max.nonone)/2
  }
  if (is.numeric(N0) == FALSE | length(N0) > 2) {
    stop("Argument N0 is not correctly specified")
  } else if (length(N0) == 1) {
    N0 = rep(N0, 2)
  }
  if (is.numeric(target.N) == FALSE) {
    stop("Argument target.N is not correctly specified")
  }
  
  ngenes = length(p.values)
  # step 2 fit mixture model using CBUM+CDD method
  Fitted.Model = MixtureModel.Fittting.pi0(p.values, thresh.p = thresh.p, restrict = F)
  
  # step3 caluclate posterior probability of being DMR
  lambda = as.numeric(Fitted.Model[1])
  r = as.numeric(Fitted.Model[2])
  s = as.numeric(Fitted.Model[3])
  posterior = lambda/(dbeta(p.values, r, s) * (1 - lambda) + lambda)  #prob of being non-DE
  
  #CDD
  lambda.cdd = convest(p.values)[[1]]
  #MLE
  lambda.mle = as.numeric(MixtureModel.Fittting.pi0.MLE(p.values, thresh.p = thresh.p, restrict = F)[[1]])
  return(list(lambda=lambda, lambda.CDD=lambda.cdd, lambda.MLE=lambda.mle, p.values=p.values, posterior=posterior))
}


####################################
# Esimtate EDR Inflate both N and R#
####################################
Estimate.EDR.from.pilot <- function(res, N0, target.N, thresh.p = 0.005, FDR = 0.05, M = 10) {
  p.values<-res$p.values
  #model<-res$model
  factorR<-res$factorR
  
  if (all(is.na(p.values))) {
    stop("all p-values were NA; nothing to compute")
  }
  if (any(is.na(p.values))) {
    p.values <- p.values[!is.na(p.values)]
  }
  if (min(p.values) == 0) {
    min.nonzero <- min(p.values[p.values > 0])
    p.values[p.values == 0] <- min.nonzero/2
  }
  if (max(p.values) == 1) {
    max.nonone <- max(p.values[p.values < 1])
    p.values[p.values == 1] <- 1 - (1 - max.nonone)/2
  }
  if (class(N0) != "numeric" | length(N0) > 2) {
    stop("Argument N0 is not correctly specified")
  } else if (length(N0) == 1) {
    N0 = rep(N0, 2)
  }
  if (!class(target.N) %in% c("matrix","numeric")) {
    stop("Argument target.N is not correctly specified")
  } else if(class(target.N)=="numeric") {
    target.N=matrix( rep( target.N, 2 ), ncol=2)
  }
  
  ngenes = length(p.values)
  # step 2 fit mixture model using CBUM+CDD method
  Fitted.Model = MixtureModel.Fittting.pi0(p.values, thresh.p = thresh.p, restrict = F)
  
  # step3 caluclate posterior probability of being DMR
  lambda = as.numeric(Fitted.Model[1])
  r = as.numeric(Fitted.Model[2])
  s = as.numeric(Fitted.Model[3])
  posterior = lambda/(dbeta(p.values, r, s) * (1 - lambda) + lambda)  #prob of being non-DE
  
  if (lambda > 0.99 | lambda < 0.5) {
    stop(paste("lambda >0.99 or <0.5,lambda=", Fitted.Model[1]))
  }
  
  parameter = list(n.old = N0, n.new = target.N)
  
  # step 4: calculate estimated EDR using parametric bootstrap method
  EDR=DeclareDMR=FDR.matrix=matrix(,nrow=ncol(factorR),ncol=nrow(target.N))
  rownames(EDR)=rownames(DeclareDMR)=rownames(FDR.matrix)=colnames(factorR)
  colnames(EDR)=colnames(DeclareDMR)=colnames(FDR.matrix)=apply(target.N, 1, function(x) paste(x, sep="", collapse=" vs "))
  for(i in 1:ncol(factorR)){
    Result = lapply(1:M, function(x) Resampling(target.N, posterior, p.values, parameter, ngenes, FDR, factorR[,i]))
    ave.result <- Reduce("+", Result)/length(Result)
    EDR[i,]<-ave.result["EDR", ]
    DeclareDMR[i,]<-round(ave.result["DeclareDMR", ])
    FDR.matrix[i,]<-ave.result["FDR", ]
  }
  pred <- list(EDR = EDR, DeclareDMR = DeclareDMR, FDR.matrix = FDR.matrix)
  class(pred) <- append(class(pred),"MethylSeqDesign")
  return(pred)
}

# calculate estimated EDR using parametric bootstrap method
Resampling <- function(target.N, posterior, p.values, parameter, ngenes, FDR, factorR) {
  
  DE_status_posterior = sapply(1:length(posterior), function(x) sample(c(FALSE, TRUE), 1, prob = c(posterior[x], 1 - posterior[x]), replace = TRUE))
  
  transform.p.values.norm <- function(p.values.each, DE_status_each, parameter, factorR.each) {
    n.old <- parameter[[1]]
    n.new <- parameter[[2]]
    
    if (DE_status_each) {
      statistic.old = qnorm(p.values.each/2, lower.tail = FALSE)
      statistic.new = statistic.old * sqrt((apply(n.new,1,prod)/prod(n.old)) * (sum(n.old)/apply(n.new,1,sum))) *sqrt(factorR.each)
      p.values.star.star = (1 - pnorm(abs(statistic.new))) * 2
      return(p.values.star.star)
    } else return(rep(p.values.each, nrow(n.new)))
  }
  
  p.values.star.posterior = sapply(1:ngenes, function(x) transform.p.values.norm(p.values[x], DE_status_posterior[x], parameter = parameter, factorR.each=factorR[x]))
  
  p.values.star.posterior = matrix(p.values.star.posterior, ncol = nrow(target.N), byrow = T)
  
  Estimate_Posterior <- function(p.values.star.star) {
    if (min(p.values.star.star, na.rm = T) == 0) {
      min.nonzero <- min(p.values.star.star[p.values.star.star > 0])
      p.values.star.star[p.values.star.star == 0] <- min.nonzero/2
    }
    if (max(p.values.star.star, na.rm = T) == 1) {
      max.non1 <- max(p.values.star.star[p.values.star.star < 1])
      p.values.star.star[p.values.star.star == 1] <- (max.non1 + 1)/2
    }
    ## empirical FDR control
    p.DE = p.values.star.star[DE_status_posterior]
    p.nonDE = p.values.star.star[!DE_status_posterior]
    p.unique = sort(unique(p.values.star.star))
    
    FDR.unique <- vector(length = length(p.unique))
    for (i in 1:length(p.unique)) {
      FDR.unique[i] = sum(p.nonDE <= p.unique[i])/sum(p.values.star.star <= p.unique[i])
    }
    index = which(FDR.unique <= FDR)
    p.values.cut = if (length(index)) {
      p.unique[max(index)]
    } else 0  #all p value cutoffs exceed the specified FDR level, fails to control FDR
    FDR.cut<-if(length(index)){
      FDR.unique[max(index)]
    } else 0
    
    if (min(p.values.star.star) > p.values.cut | is.na(p.values.cut)) {
      Declare_status = rep("nonDMR", ngenes)
    } else {
      Declare_status = DE_status_posterior
      Declare_status[which(p.values.star.star <= p.values.cut)] = "DMR"
      Declare_status[-which(p.values.star.star <= p.values.cut)] = "nonDMR"
    }
    A = sum((Declare_status == "nonDMR") * (!DE_status_posterior))
    B = sum((Declare_status == "nonDMR") * (DE_status_posterior))
    C = sum((Declare_status == "DMR") * (!DE_status_posterior))
    D = sum((Declare_status == "DMR") * (DE_status_posterior))
    if ((C + D) == 0) {
      TP_hat_post = 0
      ## no declared DMR
    } else {
      TP_hat_post = D/(C + D)
    }
    if ((A + B) == 0) {
      TN_hat_post = 0
      ## no declared DMR
    } else {
      TN_hat_post = A/(A + B)
    }
    EDR_post = D/(B + D)
    Declare_post = sum(Declare_status == "DMR")
    return(c(TP_hat_post, TN_hat_post, EDR_post, Declare_post, A = A, B = B, C = C, D = D, p.values.cut = p.values.cut, FDR = FDR.cut))
  }
  
  
  Estimate.Posterior.Result = matrix(apply(p.values.star.posterior, 2, Estimate_Posterior), ncol = nrow(target.N))
  
  Estimate.Posterior.Result <- round(as.matrix(Estimate.Posterior.Result), 3)
  
  row.names(Estimate.Posterior.Result) = c("TP", "TN", "EDR", "DeclareDMR", "A", "B", "C", "D", "FDRcut", "FDR")
  colnames(Estimate.Posterior.Result) = apply(target.N, 1, function(x) paste(x, sep="", collapse=" vs "))
  
  return(Estimate.Posterior.Result)
}

#--------------------
#  plot
#--------------------
plot2d<-function(x){
  UseMethod("plot2d", x)
}
plot3d<-function(x){
  UseMethod("plot3d", x)
}
# 2d plot
plot2d.MethylSeqDesign<-function(x){
  N<-as.numeric(sapply(strsplit(colnames(x$EDR)," vs "),function(x) x[1]))
  prop<-as.numeric(rownames(x$EDR))
  plot2d.data<-data.frame(sample.size=as.numeric(sapply(N,function(x) rep(x, length(prop)))),
                          depth.per.sample=rep(prop,length(N)),
                          estimated.EDR=as.numeric(x$EDR))
  p <- plot_ly(plot2d.data, x = ~sample.size, y = ~depth.per.sample, z = ~estimated.EDR, type="contour")
  return(p)
}

#3d plot
plot3d.MethylSeqDesign<-function(x){
  estimated.EDR<-as.matrix(x$EDR)
  p <- plot_ly(z=~estimated.EDR) %>%
    add_surface() %>% 
    layout(scene = list(xaxis = list(title = 'Sample size'),
                        yaxis = list(title = 'Seqencing depth per sample'),
                        zaxis = list(title = 'estimated EDR')))
  return(p)
}

#--------------------------------------------#
# Optimize N and R given budget or target EDR#
#--------------------------------------------#
#price_per_lane=1500
#fix_price_per_sample=300
#depth_per_lane=250
designOptim<-function(EDR.result, pilot.depth, R, N=seq(4,50,2), targetEDR= NULL, budget=NULL,
                       price_per_lane=1500, fix_price_per_sample=300, depth_per_lane=250, align_rate=0.5){
  if (!"ggplot2" %in% rownames(installed.packages())) {
    install.packages("ggplot2")
  }
  if(is.null(targetEDR) & is.null(budget)){
    stop("Either targetEDR or buget should be given")
  }
  pilot.R<-(pilot.depth/align_rate)/depth_per_lane
  R=R[pilot.R>=R]
  sample.size<-t(matrix(rep(2*N, length(R)),ncol=length(R)))
  cost.matrix<-ceiling(sweep(sample.size, 1, depth_per_lane*R, FUN="*")/depth_per_lane)*price_per_lane + sample.size*fix_price_per_sample
  # Given Power, minimize budget
  if(is.null(budget) & !is.null(targetEDR)){
    #Find the minimun budget
    min.budget<-min(cost.matrix[EDR.result$EDR>=targetEDR])
    EDR.achieve<-max(EDR.result$EDR[cost.matrix<=min.budget])
    optim.R<-R[which(cost.matrix==min.budget & EDR.result$EDR==EDR.achieve, arr.ind=TRUE)[1]]
    optim.N<-N[which(cost.matrix==min.budget & EDR.result$EDR==EDR.achieve, arr.ind=TRUE)[2]]
    optim.res=list(optim.N=optim.N, optim.R=optim.R, EDR.target=targetEDR, EDR.achive=EDR.achieve, min.budget=min.budget)
    #Find all admissible N and R
    admis.matrix<-cost.matrix
    for(i in 1:nrow(admis.matrix)){
      for(j in 1:ncol(admis.matrix)){
        admis.matrix[i,j]<-ifelse(sum(cost.matrix<cost.matrix[i,j] &  EDR.result$EDR>EDR.result$EDR[i,j])>0, 0, 1)
      }
    }
    #plot1
    cols<-ifelse(as.numeric(admis.matrix)==0, "grey", "black")
    if(0){
      plot(x=rep(N, each=length(R)), y=rep(R,length(N)), pch=4, col=cols, xlab="N", ylab="R", main="case2")
      points(x=optim.N, y=optim.R, col=2, pch=4)
      
    }
    dat<-data.frame(N=rep(N, each=length(R)),R=rep(R,length(N)))
    library(ggplot2)
    p1<-ggplot(dat, aes(x=N,y=R))+geom_point(col=cols, size=3, shape=4)+ 
      geom_point(aes(x=optim.N, y=optim.R), colour="red", size=3, shape=4)+
      ggtitle("Scenario 1")+theme_bw()+theme(plot.title = element_text(hjust = 0.5))
    
    #plot2
    cols<-ifelse(as.numeric(admis.matrix)==0, "grey", "black")
    if(0){
      plot(x=as.numeric(cost.matrix), y=as.numeric(EDR.result$EDR), pch=4, col=cols, xlab="budget", ylab="EDR", xlim=c(min.budget-10000,min.budget+10000), ylim=c(targetEDR-0.15, min(targetEDR+0.15,1)), main="case2")
      points(x=min.budget, y=EDR.achieve, pch=4, col=2)
      abline(h=targetEDR,lty=2)
    }
    dat<-data.frame(cost=as.numeric(cost.matrix),EDR=as.numeric(EDR.result$EDR))
    library(ggplot2)
    p2<-ggplot(dat, aes(x=cost,y=EDR))+geom_point(col=cols, size=3, shape=4)+
      geom_point(aes(x=min.budget, y=EDR.achieve), colour="red", size=3, shape=4)+
      geom_hline(yintercept=targetEDR,colour="blue",size=1,linetype="dashed" )+
      annotate("text", median(dat$cost), targetEDR, vjust = -1, label = paste0("target EDR = ", targetEDR))+
      #xlim(0,25000)+ylim(0.65,1)+
      ggtitle("Scenario 1")+theme_bw()+theme(plot.title = element_text(hjust = 0.5))
  }
  # Given budget, maximize power
  if(!is.null(budget) & is.null(targetEDR)){
    max.EDR<-max(EDR.result$EDR[cost.matrix<=budget])
    budget.achieve<-min(cost.matrix[EDR.result$EDR>=max.EDR])
    optim.R<-R[which(EDR.result$EDR==max.EDR & cost.matrix==budget.achieve, arr.ind=TRUE)[1]]
    optim.N<-N[which(EDR.result$EDR==max.EDR & cost.matrix==budget.achieve, arr.ind=TRUE)[2]]
    optim.res<-list(max.EDR=max.EDR, optim.N=optim.N, optim.R=optim.R)
    optim.res=list(optim.N=optim.N, optim.R=optim.R, budget=budget, cost=budget.achieve, max.EDR=max.EDR)
    #Find all admissible N and R
    admis.matrix<-cost.matrix
    for(i in 1:nrow(admis.matrix)){
      for(j in 1:ncol(admis.matrix)){
        admis.matrix[i,j]<-ifelse(sum(cost.matrix<cost.matrix[i,j] &  EDR.result$EDR>EDR.result$EDR[i,j])>0, 0, 1)
      }
    }
    #plot1
    cols<-ifelse(as.numeric(admis.matrix)==0, "grey", "black")
    if(0){
      plot(x=rep(N, each=length(R)), y=rep(R,length(N)), pch=4, col=cols, xlab="N", ylab="R", main="case1")
      points(x=optim.N, y=optim.R, col=2, pch=4)  
    }
    dat<-data.frame(N=rep(N, each=length(R)),R=rep(R,length(N)))
    library(ggplot2)
    p1<-ggplot(dat, aes(x=N,y=R))+geom_point(col=cols, size=3, shape=4)+
      geom_point(aes(x=optim.N, y=optim.R), colour="red", size=3, shape=4)+
      ggtitle("Scenario 2")+theme_bw()+theme(plot.title = element_text(hjust = 0.5))
    
    
    #plot2
    cols<-ifelse(as.numeric(admis.matrix)==0, "grey", "black")
    if(0){
      plot(x=as.numeric(cost.matrix), y=as.numeric(EDR.result$EDR), pch=4, col=cols, ylim=c(max.EDR-0.15, min((max.EDR+0.15),1)), xlim=c(budget-10000, budget+10000), xlab="budget", ylab="EDR", main="case1")
      points(x=budget.achieve, y=max.EDR, pch=4, col=2)
      abline(v=budget,lty=2)
    }
    dat<-data.frame(cost=as.numeric(cost.matrix),EDR=as.numeric(EDR.result$EDR))
    library(ggplot2)
    p2<-ggplot(dat, aes(x=cost,y=EDR))+geom_point(col=cols, size=3, shape=4)+
      geom_point(aes(x=budget.achieve, y=max.EDR), colour="red", size=3, shape=4)+
      geom_vline(xintercept=budget,colour="blue",size=1,linetype="dashed" )+
      annotate("text", budget, median(dat$EDR), vjust = -1, label = paste0("budget = ", budget))+
      #xlim(10000,30000)+ylim(0.7,1)+
      ggtitle("Scenario 2")+theme_bw()+theme(plot.title = element_text(hjust = 0.5))
  }
  return(list(res=optim.res,plot1=p1,plot2=p2))
}
