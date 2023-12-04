rm(list=ls())
setwd("C:/Users/kuan liu/Dropbox (UT_KLiu)/thesis-simulation/paper1/")
library(mvtnorm)
library(nnet)
library(MCMCpack)
library(wgeesel)
library(lme4)
library(geepack)
library(MuMIn)
library(R2jags)
library(runjags)

options(warn=-1)

expit<-function(x){1/(1+exp(-x))}
logit<-function(p){log(p)-log(1-p)}

loglik<-function(param,resp1,resp2,resp3,mmat1,mmat2,mmat3,weight){

  n<-length(resp1)
  theta <- param[1:3]

  sigma11 <- param[4]
  sigma22 <- param[5]
  sigma33 <- param[6]

  weight1<-as.vector(weight[,1])
  weight2<-as.vector(weight[,2])
  weight3<-as.vector(weight[,3])

  e1<- (resp1 - mmat1%*%theta)
  e2<- (resp2 - mmat2%*%theta)
  e3<- (resp3 - mmat3%*%theta)

  logl <- NA

  logl1<- -0.5*log(sigma11) - 0.5*((e1)^2)/(sigma11)

  logl2<- -0.5*log(sigma22) - 0.5*((e2)^2)/(sigma22)

  logl3<- -0.5*log(sigma33) - 0.5*((e3)^2)/(sigma33)

  logl<- sum(weight1*logl1)+sum(weight2*logl2)+sum(weight3*logl3)

  return(logl)

}


cat( " model {

     #N = nobs
     for (i in 1:N) {

     # conditional treatment assignment model, v2;
     z2[i] ~ dbern(p2[i])
     logit(p2[i]) <- b20 + b21*x2[i] + b22*y2[i] + b23*z1[i]

     # marginal treatment assignment model, v2;
     z2s[i] ~ dbern(p2s[i])
     logit(p2s[i]) <- bs20 + bs21*z1[i]

     # conditional treatment assignment model, v1;
     z1[i] ~ dbern(p1[i])
     logit(p1[i]) <- b10 + b11*x1[i] + b12*y1[i]

     # marginal treatment assignment model, v1;
     z1s[i] ~ dbern(p1s[i])
     logit(p1s[i]) <- bs10

     }

     # Priors
     b10 ~ dunif(-10,10) #true 0;
     b11 ~ dunif(-10,10) #true 0.2;
     b12 ~ dunif(-10,10) #true -0.05;

     b20 ~ dunif(-10,10)
     b21 ~ dunif(-10,10)
     b22 ~ dunif(-10,10)
     b23 ~ dunif(-10,10) #true 2

     bs10 ~ dunif(-10,10)
     bs20 ~ dunif(-10,10)
     bs21 ~ dunif(-10,10)


     }",
     file = "model_unif3.txt")


# Bayesian inference 1. BMSM
# first obtain MCMC sample for weights! from posterier distribution of treatment assignment parameters;
jags.data<-list(x1= obs$x1, y1=obs$y1, z1=obs$z1,  x2=obs$x2, y2=obs$y2, z2=obs$z2, z1s = obs$z1 , z2s=obs$z2, N = ntot)
jags.params<-c("b10","b11","b12","b20","b21", "b22","b23",
               "bs10","bs20","bs21")

jags.inits<-function(){list(b10=0.1,
                            b11=0.1,
                            b12=0.1,
                            b20=0.1,
                            b21=0.1,
                            b22=0.1,
                            b23 = 0.1,
                            bs10 = 0.1,
                            bs20 = 0.1,
                            bs21 = 0.1)}

jagsfit<- jags(data = list(x1= obs$x1,
                           y1=obs$y1,
                           z1=obs$z1,
                           x2=obs$x2,
                           y2=obs$y2,
                           z2=obs$z2,
                           z1s = obs$z1 ,
                           z2s=obs$z2,
                           N = ntot),
               inits = jags.inits,
               jags.params,
               n.iter = 10000,
               model.file = "model_unif3.txt",
               n.chains = 1,
               n.burnin = 5000,
               n.thin = 5)

# or to use some plots in coda
# use as.mcmmc to convert rjags object into mcmc.list
jags.mcmc <- as.mcmc(jagsfit)
out.mcmc <- as.matrix(jags.mcmc[[1]])

samplesize = 1000
ntot = 500
obs_prob1<-matrix(NA, samplesize, ntot)
exp_prob1<-matrix(NA, samplesize, ntot)
obs_prob2<-matrix(NA, samplesize, ntot)
exp_prob2<-matrix(NA, samplesize, ntot)

#calulating the MCMC weights;
for (i2 in 1:(samplesize)){
  for (j2 in 1:(ntot)){

    exp_prob1[i2,j2] <- (exp(z[j2,1]*out.mcmc[i2,8]))/(1.0+exp(out.mcmc[i2,8]))
    exp_prob2[i2,j2] <- exp_prob1[i2,j2]*(exp(z[j2,2]*(out.mcmc[i2,9]+out.mcmc[i2,10]*z[j2,1])))/(1.0+exp(out.mcmc[i2,9]+out.mcmc[i2,10]*z[j2,1]))

    obs_prob1[i2,j2] <- (exp(z[j2,1]*(out.mcmc[i2,1] + out.mcmc[i2,2]*Lobs[j2,1]+out.mcmc[i2,3]*Lobs[j2,2])))/(1.0+exp(out.mcmc[i2,1] + out.mcmc[i2,2]*Lobs[j2,1]+out.mcmc[i2,3]*Lobs[j2,2]))
    obs_prob2[i2,j2] <- obs_prob1[i2,j2]*(exp(z[j2,2]*(out.mcmc[i2,4]+out.mcmc[i2,5]*Lobs[j2,3]+out.mcmc[i2,6]*Lobs[j2,4]+out.mcmc[i2,7]*z[j2,1])))/(1.0+exp(out.mcmc[i2,4]+out.mcmc[i2,5]*Lobs[j2,3]+out.mcmc[i2,6]*Lobs[j2,4]+out.mcmc[i2,7]*z[j2,1]))

  }

  # if (i2 %% 50 == 0) {
  #   print(i2)
  # }
}

wmean2_s <- colSums(exp_prob1)/colSums(obs_prob1)
wmean3_s <- colSums(exp_prob2)/colSums(obs_prob2)
wmean_s<-cbind(rep(1,ntot),wmean2_s, wmean3_s)

#Multinomial sampling - unweighted error variance:

inits1<-c(30,5,0.2,16,36,49)
nboot = 1000

bootest1<-numeric(nboot)

for (i in 1:nboot) {
  alpha <- rep(1/ntot, ntot)
  bootidx <- as.matrix(rep(1:ntot, rmultinom(1, ntot, alpha)))

  maxim <- optim(inits1, fn=loglik,
                 resp1=obs$y1[bootidx],
                 resp2=obs$y2[bootidx],
                 resp3=obs$y3[bootidx],
                 mmat1=cbind(1,obs$cumz1[bootidx],1),
                 mmat2=cbind(1,obs$cumz2[bootidx],2),
                 mmat3=cbind(1,obs$cumz3[bootidx],3),
                 weight=wmean_s[bootidx,],
                 control=list(fnscale=-1),
                 method='BFGS', hessian=F)
  bootest1[i] <- maxim$par[2]

  # if (i %% 50 == 0) {
  #   print(i)
  # }
}

mean(bootest1)
var(bootest1)
quantile(bootest1, probs=c(0.025,0.975))

# Dirichlet sampling - unweighted error:
bootest3<-numeric(nboot)

for (j in 1:nboot) {
  alpha <- as.numeric(rdirichlet(1, rep(1.0, ntot)))

  maxim <- optim(inits1,
                 fn=loglik,
                 resp1=obs$y1,
                 resp2=obs$y2,
                 resp3=obs$y3,
                 mmat1=cbind(1,obs$cumz1,1),
                 mmat2=cbind(1,obs$cumz2,2),
                 mmat3=cbind(1,obs$cumz3,3),
                 weight=alpha*wmean_s,
                 control=list(fnscale=-1), method='BFGS', hessian=F)
  bootest3[j] <- maxim$par[2]

  # if (j %% 50 == 0) {
  #   print(j)
  # }
}


mean(bootest3)
var(bootest3)
quantile(bootest3, probs=c(0.025,0.975))
