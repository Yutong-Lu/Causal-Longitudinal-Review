rm(list=ls())

library(MCMCpack)
library(R2jags)
library(tidyverse)
options(warn=-1)

expit<-function(x){1/(1+exp(-x))}
logit<-function(p){log(p)-log(1-p)}


cat( " model {

     #N = nobs
     for (i in 1:N) {

     # conditional treatment assignment model, visit 2;
     a2[i] ~ dbern(p2[i])
     logit(p2[i]) <- b20 + b21*w1[i] + b22*w2[i] + b23*L11[i] + b24*L21[i] + b25*a1[i] + b26*L12[i] + b27*L22[i]

     # conditional treatment assignment model,visit 1;
     a1[i] ~ dbern(p1[i])
     logit(p1[i]) <- b10 + b11*w1[i] + b12*w2[i] + b13*L11[i] + b14*L21[i]

     # marginal treatment assignment model, visit 2;
     a2s[i] ~ dbern(p2s[i])
     logit(p2s[i]) <- bs20 + bs21*a1s[i]

     # marginal treatment assignment model, visit 1;
     a1s[i] ~ dbern(p1s[i])
     logit(p1s[i]) <- bs10

     }
     
     mp1s <- mean(p1s[])

     # Priors
     bs10 ~ dnorm(0,10) #intercept;
     bs20 ~ dnorm(0,10) #intercept;
     bs21 ~ dnorm(0,5) #OR;

     b10 ~ dnorm(0,10)
     b11 ~ dnorm(0,5)
     b12 ~ dnorm(0,5)
     b13 ~ dnorm(0,5)
     b14 ~ dnorm(0,5)

     b20 ~ dnorm(0,10)
     b21 ~ dnorm(0,5)
     b22 ~ dnorm(0,5)
     b23 ~ dnorm(0,5)
     b24 ~ dnorm(0,5)
     b25 ~ dnorm(0,5)
     b26 ~ dnorm(0,5)
     b27 ~ dnorm(0,5)

     }",
     file = "model_norm_testdata.txt")

# read in test analysis data;
testdata <- read_csv("testjan29_continuous_outcome_data.csv",show_col_types = FALSE)

# Assignment to Xiao, please write out the treatment assignment model given the testdata!
# this is the 1st step to make sure you know how to specify treatment assignment model;

# Bayesian inference 1. BMSM
# Step 1: Bayesian parametric estimation of treatment assignment weights;
jags.data<-list(w1=testdata$w1,
                w2=testdata$w2,
                L11=testdata$L1_1,
                L12=testdata$L1_2,
                L21=testdata$L2_1,
                L22=testdata$L2_2,
                a1 = testdata$a_1,
                a2 = testdata$a_2,
                a1s = testdata$a_1,
                a2s = testdata$a_2,
                N = length(testdata$y))

jags.params<-c("bs10","bs20","bs21",
               "b10","b11","b12","b13","b14",
               "b20","b21", "b22","b23","b24", "b25","b26","b27"
               )

jags.inits<-function(){list(bs10 = 0.1,bs20 = 0.1,bs21 = 0.1,
                            b10=0.1,b11=0.1,b12=0.1,b13=0.1,b14=0.1,
                            b20=0.1,b21=0.1,b22=0.1,b23 = 0.1,b24=0.1,b25=0.1,b26=0.1,b27 = 0.1)}

jagsfit<- jags(data = jags.data,
               inits = jags.inits,
               jags.params,
               n.iter = 10000,
               model.file = "model_norm_testdata.txt",
               n.chains = 1,
               n.burnin = 5000,
               n.thin = 5)

# or to use some plots in coda
# use as.mcmmc to convert rjags object into mcmc.list
jags.mcmc <- as.mcmc(jagsfit)
out.mcmc <- as.matrix(jags.mcmc[[1]])
geweke.diag(out.mcmc)

pdraws = (10000 - 5000)/5 # (n.iter - n.burnin)/n.thin;
n = length(testdata$y)

p1<-matrix(NA, n, pdraws)
p1s<-matrix(NA, n, pdraws)
p2<-matrix(NA, n, pdraws)
p2s<-matrix(NA, n, pdraws)

#calculating the MCMC weights;


for (j in 1:pdraws) {

# treatment assignment probability for the ith patient with jth posterior parameter draws;
p1[,j]<- expit(rowSums(as.matrix(cbind(1,testdata[,1:4]))*matrix(rep(out.mcmc[j,1:5],n),byrow = T,nrow=n)))
p2[,j]<- p1[,j]*expit(rowSums(as.matrix(cbind(1,testdata[,1:7]))*matrix(rep(out.mcmc[j,6:13],n),byrow = T,nrow=n)))

p1s[,j]<- expit(as.matrix(out.mcmc[j,14]))
p2s[,j]<- p1s[,j]*expit(rowSums(as.matrix(cbind(1,testdata[,5]))*matrix(rep(out.mcmc[j,15:16],n),byrow = T,nrow=n)))

}

wmean <- colSums(p2s)/colSums(p2)


# Step 2: Bayesian non-parametric bootstrap to calculate causal effect;
#log-likelihood for a simple linear regression where the outcome y is continous;
#we assume the residual follows a normal distribution;
#reference https://www.stat.cmu.edu/~cshalizi/mreg/15/lectures/06/lecture-06.pdf
# and https://www.ime.unicamp.br/~cnaber/optim_1.pdf
# equation 3 is the base log-likelihood function we are coding below;
# we update this log-likelihood to include the treatment weights!

# First figure out components of the likelihood;
# 1. data: dependent variable (Y, outcome) and independent variables (A, treatment variables in our case);
# 2. residuals follow a normal distribution with
# 2.1 mean: y_i - theta*a_i; theta is a parameter vector storing treatment effects;
# 2.2 variance: sigma^2, sigma is the standard deviation;
# 3. treatment weights, 1 weight per patient;

# we require wide format data, requiring 1 patient per row;
wloglik_normal<-function(param,
                        Y,
                        A,
                        weight){
  #number of observations;
  n <- length(Y)
  theta <- param[1:dim(A)[2]] #causal parameters on the mean
  #number of parameter is determined by number of treatment variables, plus intercept;
  sigma <- param[(dim(A)[2]+1)] # the remaining the parameter represent the standard deviation;
  mmat <- as.matrix(A) #design matrix of the causal outcome model, e.g., A = cbind(1, a_1, a_2);
  logl<- -0.5*log(sigma^2) - 0.5*((Y - mmat%*%theta)^2)/(sigma^2)
  wlogl<-sum(weight*logl)

  return(wlogl)
}

# we require wide format data, requiring 1 patient per row;
wloglik_binary<-function(param,
                         Y,
                         A,
                         weight){
  #number of observations;
  n <- length(Y)
  theta <- param[1:dim(A)[2]] #causal parameters on the mean
  #number of parameter is determined by number of treatment variables, plus intercept;
  mmat <- as.matrix(A) #design matrix of the causal outcome model, e.g., A = cbind(1, a_1, a_2);
  logl<- Y*(mmat%*%theta)-log(1+exp(mmat%*%theta))
  wlogl<-sum(weight*logl)
  
  return(wlogl)
}
# Dirichlet sampling

inits1<-c(0.1,0.1,0.1,0.1,4) #first few for mean parameters + 1 variance parameter
nboot <- 1000
bootest2<-numeric(nboot)

for (j in 1:nboot) {
  alpha <- as.numeric(rdirichlet(1, rep(1.0, length(testdata$y))))
  maxim <- optim(inits1,
                 fn=wloglik_binary,
                 Y=testdata$y,
                 A=cbind(1,testdata$a_1, testdata$a_2, testdata$a_1*testdata$a_2), #three mean parameters (intercept + coefficient for a_1 and coefficient for a_2);
                 weight=alpha*wmean,
                 control=list(fnscale=-1), method='BFGS', hessian=F)
  bootest2[j] <- expit(maxim$par[1]+maxim$par[2]+maxim$par[3]+maxim$par[4]) - expit(maxim$par[1]) #difference on the mean of Y between always treated and never treated;
  
  if (j %% 100 == 0) {
    print(j)
  }
}


# Dirichlet sampling

inits1<-c(0.1,0.1,0.1,0.1,4) #first few for mean parameters + 1 variance parameter
nboot <- 1000
bootest<-numeric(nboot)

for (j in 1:nboot) {
  alpha <- as.numeric(rdirichlet(1, rep(1.0, length(testdata$y))))
  maxim <- optim(inits1,
                 fn=wloglik_normal,
                 Y=testdata$y,
                 A=cbind(1,testdata$a_1, testdata$a_2, testdata$a_1*testdata$a_2), #three mean parameters (intercept + coefficient for a_1 and coefficient for a_2);
                 weight=alpha*wmean,
                 control=list(fnscale=-1), method='BFGS', hessian=F)
  bootest[j] <- maxim$par[2]+maxim$par[3]+maxim$par[4] #difference on the mean of Y between always treated and never treated;

  if (j %% 100 == 0) {
    print(j)
  }
}


mean(bootest)
# var(bootest)
sd(bootest)
quantile(bootest, probs=c(0.025,0.975))

# > mean(bootest)
# [1] -3.164371
# > var(bootest)
# [1] 0.009552388
# > quantile(bootest, probs=c(0.025,0.975))
# 2.5%     97.5%
#   -3.356955 -2.975744

#comparing to frequentist MSMs;
library(WeightIt)
Wmsm <- weightitMSM(
  list(a_1 ~ w1 + w2 + L1_1 + L2_1,
       a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
  data = testdata,
  method = "ps",
  stabilize = TRUE)

library(survey)
msm_design <- svydesign(~1, weights = Wmsm$weights, data = testdata)
fitMSM <- svyglm(y ~ a_1+a_2, design = msm_design)
summary(fitMSM)

APO_11 <- predict(fitMSM, newdata = data.frame(a_1=1,a_2=1))
APO_00 <- predict(fitMSM, newdata = data.frame(a_1=0,a_2=0))
APO_11 - APO_00
# link     SE
# 1 -3.1611 0.0758

# super similar!


