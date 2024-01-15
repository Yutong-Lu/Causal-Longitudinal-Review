###
#
###
setwd("/Users/lorrainelu/Documents/GitHub/Casual-Longitudinal-Review")
ncores<- 4
.libPaths(c(.libPaths(),"/Users/lorrainelu/Library/R/x86_64/4.1/library"))

options(warn=-1)

#parallel;
library(parallel)
library(foreach)
library(doParallel)

cl <- makeCluster(6)
registerDoParallel(cl)

mysim.cont <- function(outfile, from=1, to=4, ntot=1000, samplesize=10000) {
  
  # Simulation setup;
  # ntot=1000;
  # samplesize=10000;
  # from = 1;
  # to = 3;
  
  expit <- function(x){
    x <- exp(x)/(exp(x)+1) 
    return(x)
  }
  
  ###### BMSM ######
  
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
     a_2[i] ~ dbern(p2c[i])
     logit(p2c[i]) <- b20 + b21*w1[i]+ b22*w2[i] + b23*L1_1[i] + b24*L2_1[i] + b25*L1_2[i] + b26*L2_2[i] + b27*a_1[i]

     # marginal treatment assignment model, v2;
     a_2m[i] ~ dbern(p2m[i])
     logit(p2m[i]) <- bm20 + bm21*a_1[i]

     # conditional treatment assignment model, v1;
     a_1[i] ~ dbern(p1c[i])
     logit(p1c[i]) <- b10 + b11*w1[i]+ b12*w2[i] + b13*L1_1[i] + b14*L2_1[i]

     # marginal treatment assignment model, v1;
     a_1m[i] ~ dbern(p1m[i])
     logit(p1m[i]) <- bm10

     }

     # Priors
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
      
     bm10 ~ dnorm(0,10)
     bm20 ~ dnorm(0,10)
     bm21 ~ dnorm(0,5)


     }",
       file = "model_norm.txt")
  
  ###### BMSM ######
  
  mytrue <- function(a_1s = 1, a_2s = 1){
    #visit 1;
    w1sim <- rbinom(samplesize, 1, prob = 0.5) #50-50 male and female;
    w2sim <- rnorm(samplesize, mean = 12, sd = 4) #age;
    L1_1sim <- rbinom(samplesize, 1, prob = expit(-1 + 0.001*w1sim + 0.01*w2sim)) #a age related baseline binary clinical variable;
    L2_1sim <- rnorm(samplesize, mean = (0.001*w1sim + 0.01*w2sim), sd = 1) #a age & sex related baseline continuous clinical variable;
    
    #Visit 2, simulate all potential variables;
    L1_2sim <- rbinom(samplesize, 1, expit(-1 + 0.001*w1sim + 0.01*w2sim + 0.2*L1_1sim - 0.2*a_1s))
    L2_2sim <- rnorm(samplesize, mean = (L2_1sim - 0.001*w1sim + 0.01*w2sim - 0.2*a_1s), sd = 1)
    
    #Visit 3, simulate potential outcomes;
    truemean = (1 - 1*a_1s - 2*a_2s + 0.01*a_1s*a_2s - 0.001*w1sim + 0.01*w2sim + 0.1*L1_2sim + 0.1*L2_2sim)
    return(mean(truemean))
  }
  

  results_run<-foreach(i=from:to, .combine='rbind',.inorder=T, .verbose=T) %dopar% {
      
    .libPaths(c(.libPaths(),"/Users/lorrainelu/Library/R/x86_64/4.1/library"))
    library(tidyverse)
    library(survey)
    library(arm)
    library(ltmle)
    library(ipw) #MSM;
    library(gfoRmula) #g-computation;
    library(gtsummary)
    library(SuperLearner)
    library(WeightIt)
    library(mvtnorm)
    library(nnet)
    library(MCMCpack)
    library(wgeesel)
    library(lme4)
    library(geepack)
    library(MuMIn)
    library(R2jags)
    library(coda) #new package;
    library(runjags)
    
  # true value;
  # each method and each setting, est, se(est), low 95%CI, upper 95%CI;
    
  # data generation;
  # i=1
  set.seed(i+1234)
  
  results.it <- matrix(NA, 1, 21)
  
  # Visit 1
  w1 <- rbinom(ntot, 1, prob = 0.5) #50-50 male and female (reference);
  w2 <- rnorm(ntot, mean = 12, sd = 4) #age;
  L1_1 <- rbinom(ntot, 1, prob = expit(-1 + 0.001*w1 + 0.01*w2)) #a age & sex related baseline binary clinical variable;
  L2_1 <- rnorm(ntot, mean = (0.001*w1 + 0.01*w2), sd = 1) #a age & sex related baseline continuous clinical variable;
  
  a_1 <- rbinom(ntot, 1, expit(-1 - 0.001*w1 + 0.01*w2 - 0.1*L1_1 + 0.1*L2_1))
  
  # observational data;
  L1_2 <- rbinom(ntot, 1, expit(-1 + 0.001*w1 + 0.01*w2 + 0.2*L1_1 - 0.2*a_1))
  L2_2 <- rnorm(ntot, mean = (L2_1 - 0.001*w1 + 0.01*w2 - 0.2*a_1), sd = 1)
  a_2 <- rbinom(ntot, 1, expit(-1 - 0.1*w1 + 0.01*w2 - 0.1*L1_2 + 0.01*L2_2 + 0.5*a_1))
  
  # end-of-study outcome;
  y <- rnorm(ntot, mean = (1 - 1*a_1 - 2*a_2 + 0.01*a_1*a_2 - 0.001*w1 + 0.01*w2 + 0.1*L1_2 + 0.1*L2_2) , sd = 1)
  
  # saving final data
  dat1 <- data.frame(w1, w2, L1_1, L2_1, a_1, L1_2, L2_2, a_2, y)
  
  results.it[1,1] <- mytrue(a_1s = 1, a_2s = 1) - mytrue(a_1s = 0, a_2s = 0)

  # gformula package
  
  dat1_new <- dat1 %>%
    mutate(id = rep(1:1000)) %>% 
    pivot_longer(cols = -c(w1,w2,y,id), 
                 names_to = c("variable","visit"), 
                 names_sep = "_", 
                 values_to = "value") %>% 
    pivot_wider(names_from = variable, values_from = value) %>% 
    mutate(time = case_when(visit == 1 ~ 0,
                            visit == 2 ~ 1))
  
  dat1_new$y[dat1_new$visit == 1] <- NA
  
  id <- 'id'
  time_name <- 'time'
  # create one to for the package
  covnames <- c("L1", "L2", "a")
  outcome_name <- 'y'
  covtypes <- c('binary', 'normal', 'binary')
  histories <- c(lagged)
  histvars <- list(c('a', 'L1', 'L2'))
  
  covparams <- list(covmodels = c(L1 ~ w1 + w2 + lag1_L1 + lag1_a,
                                  L2 ~ w1 + w2 + lag1_L2 + lag1_a,
                                  a ~ w1 + w2 + lag1_L1 + L1 + lag1_L2 + L2 + lag1_a))
  ymodel <- y ~ lag1_a + a + lag1_a*a + w1 + w2 + L1 + L2 + lag1_L1 + lag1_L2
  
  intvars <- list('a', 'a')
  interventions <- list(list(c(static, rep(0, 2))),
                        list(c(static, rep(1, 2))))
  int_descript <- c('Never treat', 'Always treat')
  
  gform_cont_eof <- gformula_continuous_eof(
    obs_data = dat1_new,
    id = id,
    time_name = time_name,
    covnames =covnames,
    outcome_name = outcome_name, 
    covtypes = c("binary", "normal", "binary"),
    covparams = covparams,  ymodel = ymodel,
    intvars = intvars, interventions = interventions,
    int_descript = int_descript, 
    ref_int = 1,
    histories = c(lagged), histvars = list(c('a',"L1","L2")),
    basecovs = c("w1","w2"), 
    # boot_diag = TRUE,
    nsimul = 1000,
    nsamples = 1000, parallel = TRUE, ncores = 6,
    seed = 123)
  
  results.it[1, 2:5] <- unlist(gform_cont_eof$result[3,12:15])
  
  # ltmle packages;
  
  # ltmle without superlearner + kitchen sink gform + Qform (kitchen sink);
  tmle_model <- ltmle(dat1,
                      Anodes = c ("a_1","a_2") ,
                      Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                      Ynodes = c("y"), 
                      survivalOutcome =FALSE,
                      Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a_1",
                                 y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a_1 + a_2 + a_1*a_2"),
                      gform = c("a_1 ~ w1 + w2 + L1_1 + L2_1",
                                "a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1"),
                      abar = list(c(1,1), c(0,0)),
                      estimate.time = FALSE)
  
  out_tmle <- summary(tmle_model, estimator="tmle")
  
  results.it[1,6]<- out_tmle$effect.measures$ATE$estimate
  results.it[1,7]<- out_tmle$effect.measures$ATE$std.dev
  results.it[1,8:9] <- out_tmle$effect.measures$ATE$CI
  
  # ltmle with superlearner + kitchen sink gform + Qform;
  tmle_model_s <- ltmle(dat1,
                      Anodes = c ("a_1","a_2") ,
                      Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                      Ynodes = c("y"), 
                      survivalOutcome =FALSE,
                      Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a_1",
                                 y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a_1 + a_2 + a_1*a_2"),
                      gform = c("a_1 ~ w1 + w2 + L1_1 + L2_1",
                                "a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1"),
                      SL.library = "default", #with superlearner;
                      abar = list(c(1,1), c(0,0)),
                      estimate.time = FALSE)
  
  out_tmle_s<-summary(tmle_model_s, estimator="tmle")
  
  results.it[1,10]<- out_tmle_s$effect.measures$ATE$estimate
  results.it[1,11]<- out_tmle_s$effect.measures$ATE$std.dev
  results.it[1,12:13] <- out_tmle_s$effect.measures$ATE$CI
  
  
  # weightit package
  
  # based on simulation setup
  
  Wmsm.out.sim <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                               a_2 ~ w1 + w2 + L1_2 + L2_2 + a_1),
                          data = dat1, method = "ps",
                          stabilize = TRUE)
  
  cont_design <- svydesign(id=~1, weights = Wmsm.out.sim$weights, data = dat1)
  cont_mod <- svyglm(y ~ a_1*a_2, design = cont_design)
  p11 <- as.data.frame(predict(cont_mod, newdata = data.frame(a_1=1, a_2=1)))[[1]]
  p11.se <- as.data.frame(predict(cont_mod, newdata = data.frame(a_1=1, a_2=1)))[[2]]
  p00 <- as.data.frame(predict(cont_mod, newdata = data.frame(a_1=0, a_2=0)))[[1]]
  p00.se <- as.data.frame(predict(cont_mod, newdata = data.frame(a_1=0, a_2=0)))[[2]]
  results.it[1,14] <- p11 - p00
  results.it[1,15] <- sqrt(p11.se^2+p00.se^2)
  results.it[1,16:17] <- c(results.it[1,14]-1.96*sqrt(p11.se^2+p00.se^2),
                           results.it[1,14]+1.96*sqrt(p11.se^2+p00.se^2))
  
  
  # kitchen sink approach
  
  Wmsm.out <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                               a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                          data = dat1, method = "ps",
                          stabilize = TRUE)
  
  # estimate treatment effect
  cont_design <- svydesign(id=~1, weights = Wmsm.out$weights, data = dat1)
  cont_mod <- svyglm(y ~ a_1*a_2, design = cont_design)
  p11 <- as.data.frame(predict(cont_mod, newdata = data.frame(a_1=1, a_2=1)))[[1]]
  p11.se <- as.data.frame(predict(cont_mod, newdata = data.frame(a_1=1, a_2=1)))[[2]]
  p00 <- as.data.frame(predict(cont_mod, newdata = data.frame(a_1=0, a_2=0)))[[1]]
  p00.se <- as.data.frame(predict(cont_mod, newdata = data.frame(a_1=0, a_2=0)))[[2]]
  results.it[1,18] <- p11 - p00
  results.it[1,19] <- sqrt(p11.se^2+p00.se^2)
  results.it[1,20:21] <- c(results.it[1,18]-1.96*sqrt(p11.se^2+p00.se^2),
                           results.it[1,18]+1.96*sqrt(p11.se^2+p00.se^2))
  # Baysian MSMs
  
  # JAGS
  
  # what are the parameters in my case??
  
  # Bayesian inference 1. BMSM
  # first obtain MCMC sample for weights! from posterier distribution of treatment assignment parameters;
  jags.data<-list(w1=w1, w2=w2, a_1= a_1, a_1m= a_1, L1_1=L1_1,  L2_1=L2_1, 
                                a_2=a_2, a_2m=a_2, L1_2 = L1_2 , L2_2=L2_2, N = ntot)
  jags.params<-c("b10","b11","b12","b13","b14",
                 "b20","b21", "b22","b23","b24", "b25","b26","b27",
                 "bm10","bm20","bm21")
  
  jags.inits<-function(){list(b10=0.1,
                              b11=0.1,
                              b12=0.1,
                              b13=0.1,
                              b14=0.1,
                              b20=0.1,
                              b21=0.1,
                              b22=0.1,
                              b23 = 0.1,
                              b24=0.1,
                              b25=0.1,
                              b26 = 0.1,
                              b27 = 0.1,
                              bm10 = 0.1,
                              bm20 = 0.1,
                              bm21 = 0.1)}
  
  jagsfit<- jags(data = jags.data,
                 inits = jags.inits,
                 jags.params,
                 n.iter = 12000,
                 model.file = "model_norm.txt",
                 n.chains = 1,
                 n.burnin = 8000,
                 n.thin = 4)
  
  # or to use some plots in coda
  # use as.mcmmc to convert rjags object into mcmc.list
  jags.mcmc <- as.mcmc(jagsfit)
  out.mcmc <- as.matrix(jags.mcmc[[1]])
  
  # geweke.diag(out.mcmc)
  
  # pdraws = (n.iter-n.burnin)/n.thin
  pdraws = 1000
  # ntot = 500
  obs_prob1<-matrix(NA, pdraws, ntot)
  exp_prob1<-matrix(NA, pdraws, ntot)
  obs_prob2<-matrix(NA, pdraws, ntot)
  exp_prob2<-matrix(NA, pdraws, ntot)
  
  #calulating the MCMC weights;
  for (i2 in 1:(pdraws)){
    for (j2 in 1:(ntot)){
      
      p1m[i2,j2] <- expit(out.mcmc[j2,"bm10"])
      p1c[i2,j2]<-
        
      p2m[i2,j2] <- expit(out.mcmc[j2,"bm20"]+a_1[j2]*out.mcmc[j2,"bm21"])
      p3c[i2,j2]<-
      # exp_prob1[i2,j2] <- (exp(a_1[j2,1]*out.mcmc[i2,8]))/(1.0+exp(out.mcmc[i2,8]))
      # exp_prob2[i2,j2] <- exp_prob1[i2,j2]*(exp(z[j2,2]*(out.mcmc[i2,9]+out.mcmc[i2,10]*z[j2,1])))/(1.0+exp(out.mcmc[i2,9]+out.mcmc[i2,10]*z[j2,1]))
      
      # obs_prob1[i2,j2] <- (exp(z[j2,1]*(out.mcmc[i2,1] + out.mcmc[i2,2]*Lobs[j2,1]+out.mcmc[i2,3]*Lobs[j2,2])))/(1.0+exp(out.mcmc[i2,1] + out.mcmc[i2,2]*Lobs[j2,1]+out.mcmc[i2,3]*Lobs[j2,2]))
      # obs_prob2[i2,j2] <- obs_prob1[i2,j2]*(exp(z[j2,2]*(out.mcmc[i2,4]+out.mcmc[i2,5]*Lobs[j2,3]+out.mcmc[i2,6]*Lobs[j2,4]+out.mcmc[i2,7]*z[j2,1])))/(1.0+exp(out.mcmc[i2,4]+out.mcmc[i2,5]*Lobs[j2,3]+out.mcmc[i2,6]*Lobs[j2,4]+out.mcmc[i2,7]*z[j2,1]))
      
    }
    
    # if (i2 %% 50 == 0) {
    #   print(i2)
    # }
  }
  
  wmean2_s <- colSums(exp_prob1)/colSums(obs_prob1)
  wmean3_s <- colSums(exp_prob2)/colSums(obs_prob2)
  # wmean_s<-cbind(rep(1,ntot),wmean2_s, wmean3_s)
  
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
  
  # est.sim <- rep(NA, B)
  # est.msm <- rep(NA, B)
  # 
  # for (draw in 1:B){
  #   set.seed(draw)
  #   dat1b <- dat1[sample(1:ntot, size = ntot, replace = T),]
  #   dat1b <- tibble(dat1b)
  #   
  #   # simulation setup
  #   
  #   Wmsm.out.sim <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
  #                                    a_2 ~ w1 + w2 + L1_2 + L2_2 + a_1),
  #                               data = dat1, method = "ps",
  #                               stabilize = TRUE)
  #   
  #   cont_design <- svydesign(id=~1, weights = Wmsm.out.sim$weights, data = dat1b)
  #   cont_mod <- svyglm(y ~ a_1*a_2, design = cont_design)
  #   p11 <- predict(cont_mod, newdata = data.frame(a_1=1, a_2=1))[1]
  #   p00 <- predict(cont_mod, newdata = data.frame(a_1=0, a_2=0))[1]
  #   est.sim[draw] <- p11-p00
  #   
  #   # kitchen sink approach
  #   
  #   Wmsm.out <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
  #                                a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
  #                           data = dat1, method = "ps",
  #                           stabilize = TRUE)
  #   
  #   # estimate treatment effect
  #   cont_design <- svydesign(id=~1, weights = Wmsm.out$weights, data = dat1b)
  #   cont_mod <- svyglm(y ~ a_1*a_2, design = cont_design)
  #   p11 <- predict(cont_mod, newdata = data.frame(a_1=1, a_2=1))[1]
  #   p00 <- predict(cont_mod, newdata = data.frame(a_1=0, a_2=0))[1]
  #   est.msm[draw] <- p11-p00
  #   
  # }
  # 
  # est.sim <- unlist(est.sim)
  # # results.it[1,14]<-mean(est.sim)
  # results.it[1,15]<-sd(est.sim)
  # results.it[1,16:17]<-c(results.it[1,14]-1.96*sd(est.sim), results.it[1,14]+1.96*sd(est.sim))
  # 
  # est.msm <- unlist(est.msm) # unlist to calculate mean and sd
  # # results.it[1,18]<-mean(est.msm)
  # results.it[1,19]<-sd(est.msm)
  # results.it[1,20:21]<-c(results.it[1,18]-1.96*sd(est.msm), results.it[1,18]+1.96*sd(est.msm))
  
  cbind(i,results.it)
  
}

# outfile <-"textcontsim"
  
write.table(results_run, file = paste0(outfile,".txt"), row.names = FALSE, col.names = FALSE)

}

start_time <-Sys.time()
mysim.cont(paste0("Apr29","cont","run1"), from=1, to=1000, ntot=1000, samplesize=10000)
end_time <- Sys.time()
end_time - start_time


