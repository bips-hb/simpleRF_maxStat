#### Problems.R #####
#### functions to generate simulated data

#problem requires parameters data and job

#1) Simple Simulation ####
SimSimple <-function(data, job, nobs=100, p=10, confmean=2, confsd=1, beta=1, sd=1){
  #Null case: single continuous confounder z1 that affects outcome and the first predictor x1

  #*simulate confounders
  #simulate normally distributed confounder z1
  z1 <- rnorm(nobs, mean=confmean, sd=confsd)
  #*simulate predictors
  #predictor x1 (affected by confounder z1)
  x1 <- rnorm(n=nobs, mean=beta*z1, sd=sd)
  #unconfounded predictors x2-x10 (assumed to be independent & same variance as x1)
  X <- data.frame(sapply(1:(p-1), function(x) rnorm(nobs,mean=0,sd=sqrt(beta^2*confsd^2+sd^2))))
  colnames(X)<-paste('x', 2:p, sep='')
  #*simulate outcome (affected by confounder z1)
  y <- rnorm(n=nobs, mean=beta*z1, sd=sd)
  #return data set
  return(data.frame(cbind(y=y, z1=z1, x1=x1, X)))
}

SimSimplePower <-function(data, job, nobs=100, p=10, confmean=2, confsd=1, beta=1, sd=1, betax2=0.5){
  #Power case: single continuous confounder z1 that affects outcome and the first predictor x1
  #& additionally outcome affected by x2
  
  #*simulate confounders
  #simulate normally distributed confounder z1
  z1 <- rnorm(nobs, mean=confmean, sd=confsd)
  #*simulate predictors
  #predictor x1 (affected by confounder z1)
  x1 <- rnorm(n=nobs, mean=beta*z1, sd=sd)
  #unconfounded predictors x2-x10 (assumed to be independent & same variance as x1)
  X <- data.frame(sapply(1:(p-1), function(x) rnorm(nobs,mean=0,sd=sqrt(beta^2*confsd^2+sd^2))))
  colnames(X)<-paste('x', 2:p, sep='')
  #*simulate outcome (affected by confounder z1 & predictor x2)
  y <- rnorm(n=nobs, mean=beta*z1+betax2*X$x2, sd=sd)
  #return data set
  return(data.frame(cbind(y=y, z1=z1, x1=x1, X)))
}

SimSimpleLocal <- function(data, job, nobs=100, p=10, confmean=2, confsd=1, beta=1, sd=1, addx2=FALSE){
  #Local simulation: different confounder effect for subgroups defined by binary splitting variable x2
  
  #*simulate confounders
  #simulate normally distributed confounder z1
  z1 <- rnorm(nobs, mean=confmean, sd=confsd)
  #*simulate predictors
  #subgroups defined by x2
  x2 <- rbinom(nobs, size=1, prob=0.5)
  #unconfounded predictors x3-x10 (assumed to be independent & same variance as x1)
  X <- data.frame(sapply(1:(p-2), function(x) rnorm(nobs,mean=0,sd=sqrt(beta^2*confsd^2+sd^2))))
  colnames(X)<-paste('x',3:p, sep='')
  #confounded predictor x1 (affected by confounder z1 - different effect for subgroups defined by x2)
  sign <- rep(1, times=nobs)
  sign[which(x2==0)] <- -1
  x1 <- rnorm(n=nobs, mean=sign*beta*z1, sd=sd)
  #*simulate outcome (affected by z1 - different effect for subgroups defined by x2)
  #if addx2=TRUE x2 has an additional linear effect on the outcome
  if (addx2==FALSE){
    y <- rnorm(n=nobs, mean=sign*beta*z1, sd=sd)
  } else {
    y <- rnorm(n=nobs, mean=sign*beta*z1+beta*x2, sd=sd)
  }
  #return data set
  return(data.frame(cbind(y=y, z1=z1, x1=x1, x2=x2, X)))
}

#*Training & Test data ----
#create training and test data each with nobs observations for the different simulation scenarios
Pred_SimSimple <- function(data, job, nobs=100, p=10, confmean=2, confsd=1, beta=1, sd=1){
  source("Problems.R")
  #Training data
  train <- SimSimple(data, job, nobs=nobs, p=p, confmean=confmean, confsd=confsd, beta=beta, sd=sd)
  #Test data
  test <- SimSimple(data, job, nobs=nobs, p=p, confmean=confmean, confsd=confsd, beta=beta, sd=sd)
  #return training and test data
  return(list(train=train, test=test))
}

Pred_SimSimplePower <- function(data, job, nobs=100, p=10, confmean=2, confsd=1, beta=1, sd=1, betax2=0.5){
  source("Problems.R")
  #Training data
  train <- SimSimplePower(data, job, nobs=nobs, p=p, confmean=confmean, confsd=confsd, beta=beta, sd=sd, betax2=betax2)
  #Test data
  test <- SimSimplePower(data, job, nobs=nobs, p=p, confmean=confmean, confsd=confsd, beta=beta, sd=sd, betax2=betax2)
  #return training and test data
  return(list(train=train, test=test))
}

Pred_SimSimpleLocal <-function(data, job, nobs=100, p=10, confmean=2, confsd=1, beta=1, sd=1, addx2=FALSE){
  source("Problems.R")
  #Training data
  train <- SimSimpleLocal(data, job, nobs=nobs, p=p, confmean=confmean, confsd=confsd, beta=beta, sd=sd, addx2=addx2)
  #Test data
  test <- SimSimpleLocal(data, job, nobs=nobs, p=p, confmean=confmean, confsd=confsd, beta=beta, sd=sd, addx2=addx2)
  #return training and test data
  return(list(train=train, test=test))
}


#2) Genetic Simulation ####
SimSNP <- function(data, job, nobs=100, p=100, beta=1, beta0=-57, maf=0.25){
  #Null case: continuous confounder z1 & binary confounder z2 that 
  #affect outcome and the first predictor X1
  
  #*simulate confounders
  #simulate normally distributed confounder z1
  z1 <- rnorm(nobs, mean=50, sd=10)
  #simulate binary confounder z2 (bernoulli-distributed)
  z2 <- rbinom(nobs, size=1, prob=0.5)
  #*simulate predictors according to Hardy-Weinberg equilibrium:
  #prob(SNP=0)=(1-maf)^2, prob(SNP=1)=2*(1-maf)*maf, prob(SNP=2)=maf^2
  #predictors x1-x5 (affected by confounder z1 & z2)
  pi <- exp(beta0+beta*z1+beta*z2)/(1+exp(beta0+beta*z1+beta*z2)) #this is the same as plogis(linpred)
  x1 <- sapply(pi, function(x) sample(c(0,1,2), size=1, replace=TRUE, 
                                      prob=c((1-x)^2, 2*(1-x)*x, x^2)))
  x2 <- sapply(pi, function(x) sample(c(0,1,2), size=1, replace=TRUE, 
                                                prob=c((1-x)^2, 2*(1-x)*x, x^2)))
  x3 <- sapply(pi, function(x) sample(c(0,1,2), size=1, replace=TRUE, 
                                          prob=c((1-x)^2, 2*(1-x)*x, x^2)))
  x4 <- sapply(pi, function(x) sample(c(0,1,2), size=1, replace=TRUE, 
                                                prob=c((1-x)^2, 2*(1-x)*x, x^2)))
  x5 <- sapply(pi, function(x) sample(c(0,1,2), size=1, replace=TRUE, 
                                                prob=c((1-x)^2, 2*(1-x)*x, x^2)))
  #unconfounded predictors x6-x100 (assumed to be independent)
  X <- data.frame(sapply(1:(p-5), function(x) sample(c(0,1,2), nobs, replace=TRUE, 
                     prob=c((1-maf)^2, 2*(1-maf)*maf, maf^2)))) 
  colnames(X)<-paste('x', 6:p, sep='')
  #*simulate outcome (affected by confounders z1 & z2)
  y <- rnorm(n=nobs, mean=beta*z1+beta*z2, sd=10)
  #return data set
  return(data.frame(cbind(y=y, z1=z1, z2=z2, x1=x1, x2=x2, x3=x3, x4=x4, x5=x5, X)))
}

SimSNPPower <- function(data, job, nobs=100, p=100, beta=1, beta0=-57, maf=0.25, betax2=20){
  #Power case: continuous confounder z1 & binary confounder z2 that 
  #affect outcome and the first predictor x1 & additionally outcome affected by x6-x10
  
  #*simulate confounders
  #simulate normally distributed confounder z1
  z1 <- rnorm(nobs, mean=50, sd=10)
  #simulate binary confounder z2 (bernoulli-distributed)
  z2 <- rbinom(nobs, size=1, prob=0.5)
  #*simulate predictors according to Hardy-Weinberg equilibrium:
  #prob(SNP=0)=(1-maf)^2, prob(SNP=1)=2*(1-maf)*maf, prob(SNP=2)=maf^2
  #predictor x1-x5 (affected by confounder z1 & z2)
  pi <- exp(beta0+beta*z1+beta*z2)/(1+exp(beta0+beta*z1+beta*z2)) #this is the same as plogis(linpred)
  x1 <- sapply(pi, function(x) sample(c(0,1,2), size=1, replace=TRUE, 
                                                prob=c((1-x)^2, 2*(1-x)*x, x^2)))
  x2 <- sapply(pi, function(x) sample(c(0,1,2), size=1, replace=TRUE, 
                                                prob=c((1-x)^2, 2*(1-x)*x, x^2)))
  x3 <- sapply(pi, function(x) sample(c(0,1,2), size=1, replace=TRUE, 
                                                prob=c((1-x)^2, 2*(1-x)*x, x^2)))
  x4 <- sapply(pi, function(x) sample(c(0,1,2), size=1, replace=TRUE, 
                                                prob=c((1-x)^2, 2*(1-x)*x, x^2)))
  x5 <- sapply(pi, function(x) sample(c(0,1,2), size=1, replace=TRUE,                                             
                                                prob=c((1-x)^2, 2*(1-x)*x, x^2)))
  #unconfounded predictors x6-x100 (assumed to be independent)
  X <- data.frame(sapply(1:(p-5), function(x) sample(c(0,1,2), nobs, replace=TRUE, 
                                                     prob=c((1-maf)^2, 2*(1-maf)*maf, maf^2)))) 
  colnames(X)<-paste('x', 6:p, sep='')
  #*simulate outcome (affected by confounders z1 & z2 & predictors x6-x10)
  y <- rnorm(n=nobs, mean=beta*z1+beta*z2+betax2*X$x6+betax2*X$x7+betax2*X$x8+betax2*X$x9+betax2*X$x10, sd=10)
  #return data set
  return(data.frame(cbind(y=y, z1=z1, z2=z2, x1=x1, x2=x2, x3=x3, x4=x4, x5=x5, X)))
}

#*Training & Test data ----
#create training and test data each with nobs observations for the different simulation scenarios
Pred_SimSNP <-function(data, job, nobs=100, p=100, beta=1, beta0=-57, maf=0.25){
  source("Problems.R")
  #Training data
  train <- SimSNP(data, job, nobs=nobs, p=p, beta=beta, beta0=beta0, maf=maf)
  #Test data
  test <- SimSNP(data, job, nobs=nobs, p=p, beta=beta, beta0=beta0, maf=maf)
  #return training and test data
  return(list(train=train, test=test))
}

Pred_SimSNPPower <-function(data, job, nobs=100, p=100, beta=1, beta0=-57, maf=0.25, betax2=20){
  source("Problems.R")
  #Training data
  train <- SimSNPPower(data, job, nobs=nobs, p=p, beta=beta, beta0=beta0, maf=maf, betax2=betax2)
  #Test data
  test <- SimSNPPower(data, job, nobs=nobs, p=p, beta=beta, beta0=beta0, maf=maf, betax2=betax2)
  #return training and test data
  return(list(train=train, test=test))
}




