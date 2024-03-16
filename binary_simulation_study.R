###
#
###
setwd("/Users/lorrainelu/Documents/GitHub/Casual-Longitudinal-Review")
.libPaths(c(.libPaths(),"/Users/lorrainelu/Library/R/x86_64/4.1/library"))

options(warn=-1)

#parallel;
library(parallel)
library(foreach)
library(doParallel)

cl <- makeCluster(6)
registerDoParallel(cl)

mysim.bin <- function(outfile, from=1, to=4, ntot=1000, samplesize=10000, B=1000) {
  
  # Simulation setup;
  # ntot=1000;
  # samplesize=10000;
  # from = 1;
  # to = 3;
  
  library(doParallel)
  # registerDoParallel(ncores)
    
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
    # calc odds of y=1
    true_prob = expit(- 1 + log(0.9)*a_1s + log(0.85)*a_2s + log(1.01)*a_1s*a_2s
                      + log(0.95)*w1sim + log(1.05)*w2sim + log(1.15)*L1_2sim + log(1.1)*L2_2sim)
    y <- (runif(samplesize) < true_prob)
    return(y)
  }

  results_run<-foreach(i=from:to, .combine='rbind',.inorder=T, .verbose=T) %dopar% {
    
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
    library(wgeesel)
    library(lme4)
    library(geepack)
    library(MuMIn)
    library(MCMCpack)
    library(R2jags)
    library(coda) #new package;
    library(runjags)
  
    
    
    # data generation;
    set.seed(i+123)
    
    results.it <- matrix(NA, 1, 25)
  
    # Visit 1;
    w1 <- rbinom(ntot, 1, prob = 0.5) #50-50 male and female (reference);
    w2 <- rnorm(ntot, mean = 12, sd = 4) #age;
    L1_1 <- rbinom(ntot, 1, prob = expit(-1 + 0.001*w1 + 0.01*w2)) #a age & sex related baseline binary clinical variable;
    L2_1 <- rnorm(ntot, mean = (0.001*w1 + 0.01*w2), sd = 1) #a age & sex related baseline continuous clinical variable;
    
    a_1 <- rbinom(ntot, 1, expit(-1 - 0.1*w1 + 0.01*w2 - 0.1*L1_1 + 0.01*L2_1)); #exposure;
    
    # observational data;
    L1_2 <- rbinom(ntot, 1, expit(-1 + 0.001*w1 + 0.01*w2 + 0.2*L1_1 - 0.2*a_1))
    L2_2 <- rnorm(ntot, mean = (L2_1 - 0.001*w1 + 0.01*w2 - 0.2*a_1), sd = 1)
    a_2 <- rbinom(ntot, 1, expit(-1 - 0.1*w1 + 0.01*w2 - 0.1*L1_2 + 0.01*L2_2 + 0.5*a_1))
    # table(a_1, a_2)
    
    # end-of-study outcome;
    y <- rbinom(ntot, 1, prob = expit(- 1 + log(0.9)*a_1 + log(0.85)*a_2 + log(1.01)*a_1*a_2 
                                      + log(0.95)*w1 + log(1.05)*w2 + log(1.15)*L1_2 + log(1.1)*L2_2))
    
    # saving final data;
    dat2 <- data.frame(w1, w2, L1_1, L2_1, a_1, L1_2, L2_2, a_2, y)
    
    
    # getting true OR from glm model using psuedo-population;
    ylong <- c(mytrue(a_1s = 0, a_2s = 0),mytrue(a_1s = 0, a_2s = 1),mytrue(a_1s = 1, a_2s = 0),mytrue(a_1s = 1, a_2s = 1))
    z1long <- c(rep(0, 2*samplesize), rep(1,2*samplesize))
    z2long <- c(rep(0, samplesize), rep(1,samplesize), rep(0, samplesize), rep(1, samplesize))
    
    truemodel <- glm(ylong ~ z1long*z2long)
    # summary(truemodel)
    prob_11 <-predict(truemodel, newdata = data.frame(z1long=1, z2long=1), type="response")
    prob_00 <-predict(truemodel, newdata = data.frame(z1long=0, z2long=0), type="response")
    
    results.it[1,1]<-(prob_11/(1-prob_11))/(prob_00/(1-prob_00))
      
    # library(cobalt) #package to assess covariates balance by treatment;
    # 
    # #covariates balance at each visit;
    # bal.tab(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
    #              a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
    #         data = dat2, 
    #         int = FALSE,
    #         poly = 1, 
    #         estimand = "ATE", 
    #         stats = c("m"),
    #         thresholds = c(m = 0.1),
    #         which.time = .all)
    
    
    # convert prob into odds
    # N0TE: plogis and expit are very close but not exactly the same
    # prob_11<-mytrue(a_1s = 1, a_2s = 1) #potential outcome Y(a_1=1, a_2=1);
    # prob_00<-mytrue(a_1s = 0, a_2s = 0)
    
    # results.it[1,1]<-mean((prob_11/(1-prob_11))/(prob_00/(1-prob_00)))
    # results.it[1,2]<-(mean(prob_11)/(1-mean(prob_11)))/(mean(prob_00)/(1-mean(prob_00)))
    # results.it[1,3]<-mean((prob_11/(1-prob_11)))/mean((prob_00/(1-prob_00)))
    # gfoRmula package;
    # creating long format data;
    
    dat2_new <- dat2 %>%
      mutate(row = row_number()) %>%
      pivot_longer(cols = -c(w1,w2,y,row),
                   names_to = c("variable","visit"),
                   names_sep = "_",
                   values_to = "value") %>%
      pivot_wider(names_from = variable, values_from = value) %>%
      mutate(time = case_when(visit == 1 ~ 0,
                              visit == 2 ~ 1)) # time has to start with 0
  
    dat2_new$y[dat2_new$visit == 1] <- NA
  
    id <- 'row' # changed from id
    time_name <- 'time'
    covnames <- c("L1", "L2", "a")
    outcome_name <- 'y'
    covtypes <- c('binary', 'normal', 'binary')
    histories <- c(lagged, cumavg)
    histvars <- list(c('a', 'L1', 'L2'), c('L1', 'L2'))
  
    # have to specify each model for time dep cov, trt and outcome
    # need to parametrically specify all the time dep component
    covparams <- list(covmodels = c(L1 ~ w1 + w2 + lag1_L1 + lag1_a,
                                    L2 ~ w1 + w2 + lag1_L2 + lag1_a,
                                    a ~ w1 + w2 + lag1_L1 + L1 + lag1_L2 + L2 + lag1_a))
    ymodel <- y ~ lag1_a*a + w1 + w2 + L1 + L2 + lag1_L1 + lag1_L2
  
    intvars <- list('a', 'a')
  
    interventions <- list(list(c(static, rep(0, 2))),
                          list(c(static, rep(1, 2))))
    int_descript <- c('Never treat', 'Always treat')
  
    gform_bin_eof <- gformula_binary_eof(
      obs_data = dat2_new,
      id = id,
      time_name = time_name,
      covnames = covnames,
      outcome_name = outcome_name,
      covtypes = c("binary", "normal", "binary"),
      covparams = covparams,
      ymodel = ymodel,
      intvars = intvars, interventions = interventions,
      int_descript = int_descript, ref_int = 1,
      histories = c(lagged), histvars = list(c('a',"L1","L2")),
      basecovs = c("w1","w2"),
      seed=123)
  
    est.prob_00<-unlist(summary(gform_bin_eof)$result[2,4]) # need to unlist here
    est.prob_11<-unlist(summary(gform_bin_eof)$result[3,4]) # need to unlist here
    results.it[1,2] <- (est.prob_11/(1-est.prob_11))/(est.prob_00/(1-est.prob_00))
  
    est.gcomp.or <- rep(NA,B)
  
    for (draw in 1:B){
      set.seed(draw)
      dat2b <- dat2[sample(1:ntot, size = ntot, replace = T),]
      # dat2b <- tibble(dat2b)
  
      # creating long format data;
      dat2_new_b <- dat2b %>%
        mutate(row = row_number()) %>%
        pivot_longer(cols = -c(w1,w2,y,row),
                     names_to = c("variable","visit"),
                     names_sep = "_",
                     values_to = "value") %>%
        pivot_wider(names_from = variable, values_from = value) %>%
        mutate(time = case_when(visit == 1 ~ 0,
                                visit == 2 ~ 1)) # time has to start with 0
  
      dat2_new_b$y[dat2_new_b$visit == 1] <- NA
  
      id <- 'row' # changed from id
      time_name <- 'time'
      covnames <- c("L1", "L2", "a")
      outcome_name <- 'y'
      covtypes <- c('binary', 'normal', 'binary')
      histories <- c(lagged, cumavg)
      histvars <- list(c('a', 'L1', 'L2'), c('L1', 'L2'))
  
      # have to specify each model for time dep cov, trt and outcome
      # need to parametrically specify all the time dep component
      covparams <- list(covmodels = c(L1 ~ w1 + w2 + lag1_L1 + lag1_a,
                                      L2 ~ w1 + w2 + lag1_L2 + lag1_a,
                                      a ~ w1 + w2 + lag1_L1 + L1 + lag1_L2 + L2 + lag1_a))
      ymodel <- y ~ lag1_a*a + w1 + w2 + L1 + L2 + lag1_L1 + lag1_L2
  
      intvars <- list('a', 'a')
  
      interventions <- list(list(c(static, rep(0, 2))),
                            list(c(static, rep(1, 2))))
      int_descript <- c('Never treat', 'Always treat')
  
      gform_bin_eof <- gformula_binary_eof(
        obs_data = dat2_new_b,
        id = id,
        time_name = time_name,
        covnames = covnames,
        outcome_name = outcome_name,
        covtypes = c("binary", "normal", "binary"),
        covparams = covparams,
        ymodel = ymodel,
        intvars = intvars, interventions = interventions,
        int_descript = int_descript, ref_int = 1,
        histories = c(lagged), histvars = list(c('a',"L1","L2")),
        basecovs = c("w1","w2"),
        seed=123)
  
      est.prob_00<- summary(gform_bin_eof)$result[2,4]
      est.prob_11<- summary(gform_bin_eof)$result[3,4]
      est.gcomp.or[draw] <- (est.prob_11/(1-est.prob_11))/(est.prob_00/(1-est.prob_00))
  
    }
  
    est.gcomp.or <- unlist(est.gcomp.or) # unlist to calculate mean and sd
  
    # results.it[1,2] <- mean(est.gcomp.or); # within the CI OR+1.96*OR_SE;
    results.it[1,3] <- sd(est.gcomp.or)
    results.it[1,4:5] <- c(results.it[1,2]-1.96*sd(est.gcomp.or), results.it[1,2]+1.96*sd(est.gcomp.or)) # to be corrected;
  
    # ltmle package;
  
    tmle_model_noSL <- ltmle(dat2,
                             Anodes = c ("a_1","a_2") ,
                             Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"),
                             Ynodes = c("y"),
                             survivalOutcome =FALSE,
                             Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a_1",
                                        y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a_1 + a_2 + a_1*a_2"),
                             gform = c("a_1 ~ w1 + w2 + L1_1 + L2_1",
                                       "a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1"),
                             # gcomp = TRUE,
                             # iptw.only = FALSE,
                             # variance.method = "tmle",
                             # SL.library = "default", #with superlearner;
                             abar = list(c(1,1), c(0,0)),
                             estimate.time = FALSE)
  
    out_tmle <- summary(tmle_model_noSL, estimator="tmle")
  
    results.it[1,6]<- out_tmle$effect.measures$OR$estimate
    results.it[1,7]<- out_tmle$effect.measures$OR$std.dev
    results.it[1,8:9] <- out_tmle$effect.measures$OR$CI
  
    # tmle with superlearner + kitchen sink gform + Qform;
    tmle_model_SL <- ltmle(dat2,
                           Anodes = c ("a_1","a_2") ,
                           Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"),
                           Ynodes = c("y"),
                           survivalOutcome =FALSE,
                           Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a_1",
                                      y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a_1 + a_2 + a_1*a_2"),
                           gform = c("a_1 ~ w1 + w2 + L1_1 + L2_1",
                                     "a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1"),
                           # gcomp = TRUE,
                           # iptw.only = FALSE,
                           # variance.method = "tmle",
                           SL.library = "default",
                           abar = list(c(1,1), c(0,0)),
                           estimate.time = FALSE)
  
    out_tmle_s<-summary(tmle_model_SL, estimator="tmle")
  
    results.it[1,10]<- out_tmle_s$effect.measures$OR$estimate
    results.it[1,11]<- out_tmle_s$effect.measures$OR$std.dev
    results.it[1,12:13] <- out_tmle_s$effect.measures$OR$CI
  
    # weightit package
    # mean from package function, variance from bootstrap;
  
    # simulation setup
  
    Wmsm.out.sim <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                     a_2 ~ w1 + w2 + L1_2 + L2_2 + a_1),
                                data = dat2, method = "ps",
                                stabilize = TRUE)
  
    bin_design <- svydesign(id=~1, weights = Wmsm.out.sim$weights, data = dat2)
    bin_mod <- svyglm(y ~ a_1*a_2, family = "binomial", design = bin_design)
    p11 <- predict(bin_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
    p00 <- predict(bin_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
    results.it[1,14] <- (p11/(1-p11))/(p00/(1-p00))
  
    # kitchen sink
  
    Wmsm.out.sim <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                     a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                                data = dat2, method = "ps",
                                stabilize = TRUE)
  
    bin_design <- svydesign(id=~1, weights = Wmsm.out.sim$weights, data = dat2)
    bin_mod <- svyglm(y ~ a_1*a_2, family = "binomial", design = bin_design)
    p11 <- predict(bin_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
    p00 <- predict(bin_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
    results.it[1,18] <- (p11/(1-p11))/(p00/(1-p00))
    
    # try cbps and bart methods
    
    # cbps
    
    # based on simulation setup
    
    Wmsm.out.sim.cbps <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                          a_2 ~ w1 + w2 + L1_2 + L2_2 + a_1),
                                     data = dat1, method = "cbps",
                                     stabilize = TRUE)
    
    bin_design <- svydesign(id=~1, weights = Wmsm.out.sim.cbps$weights, data = dat2)
    bin_mod <- svyglm(y ~ a_1*a_2, family = "binomial", design = bin_design)
    p11 <- predict(bin_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
    p00 <- predict(bin_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
    results.it[1,22] <- (p11/(1-p11))/(p00/(1-p00))
    
    # kitchen sink approach
    
    Wmsm.out.ks.cbps <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                      a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                                 data = dat1, method = "cbps",
                                 stabilize = TRUE)
    
    # estimate treatment effect
    
    bin_design <- svydesign(id=~1, weights = Wmsm.out.ks.cbps$weights, data = dat2)
    bin_mod <- svyglm(y ~ a_1*a_2, family = "binomial", design = bin_design)
    p11 <- predict(bin_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
    p00 <- predict(bin_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
    results.it[1,26] <- (p11/(1-p11))/(p00/(1-p00))
    
    # bart
    
    # based on simulation setup
    
    Wmsm.out.sim.bart <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                          a_2 ~ w1 + w2 + L1_2 + L2_2 + a_1),
                                     data = dat1, method = "bart",
                                     stabilize = TRUE)
    
    bin_design <- svydesign(id=~1, weights = Wmsm.out.sim.bart$weights, data = dat2)
    bin_mod <- svyglm(y ~ a_1*a_2, family = "binomial", design = bin_design)
    p11 <- predict(bin_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
    p00 <- predict(bin_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
    results.it[1,30] <- (p11/(1-p11))/(p00/(1-p00))
    
    # kitchen sink approach
    
    Wmsm.out.ks.bart <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                      a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                                 data = dat1, method = "bart",
                                 stabilize = TRUE)
    
    # estimate treatment effect
    bin_design <- svydesign(id=~1, weights = Wmsm.out.ks.bart$weights, data = dat2)
    bin_mod <- svyglm(y ~ a_1*a_2, family = "binomial", design = bin_design)
    p11 <- predict(bin_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
    p00 <- predict(bin_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
    results.it[1,34] <- (p11/(1-p11))/(p00/(1-p00))
  
    est.ps.sim.or <- rep(NA, B)
    est.ps.ks.or <- rep(NA, B)
    est.cbps.sim.or <- rep(NA, B)
    est.cbps.ks.or <- rep(NA, B)
    est.bart.sim.or <- rep(NA, B)
    est.bart.ks.or <- rep(NA, B)
  
    for (draw in 1:B){
      set.seed(draw)
      dat2b <- dat2[sample(1:ntot, size = ntot, replace = T),]
      dat2b <- tibble(dat2b)
  
      # dat2b <- dat2b %>%
      #   mutate(id = rep(1:1000),
      #          cum_a = a_1 + a_2)
      
      # method ps
  
      # calculate using simulation setup
    
      Wmsm.out.sim <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                       a_2 ~ w1 + w2 + L1_2 + L2_2 + a_1),
                                  data = dat2b, method = "ps",
                                  stabilize = TRUE)
    
      bin_design <- svydesign(id=~1, weights = Wmsm.out.sim$weights, data = dat2b)
      bin_mod <- svyglm(y ~ a_1*a_2, family = "binomial", design = bin_design)
      p11 <- predict(bin_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
      p00 <- predict(bin_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
      est.ps.sim.or[draw] <- (p11/(1-p11))/(p00/(1-p00))
    
      # calculate using kitchen sink approach
    
      Wmsm.out <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                   a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                              data = dat2b, method = "ps",
                              stabilize = TRUE)
    
      d.w.msm <- svydesign(~1, weights = Wmsm.out$weights, data = dat2b)
      cum.fit <- svyglm(y ~ a_1*a_2, design = d.w.msm, family = "binomial")
      p11 <- predict(cum.fit, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
      p00 <- predict(cum.fit, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
      est.ps.ks.or[draw] <- (p11/(1-p11))/(p00/(1-p00))
      
      # method cbps
      
      Wmsm.out.sim.cbps <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                            a_2 ~ w1 + w2 + L1_2 + L2_2 + a_1),
                                       data = dat1, method = "cbps",
                                       stabilize = TRUE)
      
      bin_design <- svydesign(id=~1, weights = Wmsm.out.sim.cbps$weights, data = dat2b)
      bin_mod <- svyglm(y ~ a_1*a_2,  family = "binomial", design = bin_design)
      p11 <- predict(bin_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
      p00 <- predict(bin_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
      est.cbps.sim.or[draw] <- (p11/(1-p11))/(p00/(1-p00))
      
      # kitchen sink approach
      
      Wmsm.out.ks.cbps <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                        a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                                   data = dat1, method = "cbps",
                                   stabilize = TRUE)
      
      # estimate treatment effect
      
      bin_design <- svydesign(id=~1, weights = Wmsm.out.ks.cbps$weights, data = dat2b)
      bin_mod <- svyglm(y ~ a_1*a_2,  family = "binomial", design = bin_design)
      p11 <- predict(bin_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
      p00 <- predict(bin_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
      est.cbps.ks.or[draw] <- (p11/(1-p11))/(p00/(1-p00))
      
      # method bart
      
      Wmsm.out.sim.bart <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                            a_2 ~ w1 + w2 + L1_2 + L2_2 + a_1),
                                       data = dat1, method = "cbps",
                                       stabilize = TRUE)
      
      bin_design <- svydesign(id=~1, weights = Wmsm.out.sim.bart$weights, data = dat2b)
      bin_mod <- svyglm(y ~ a_1*a_2,  family = "binomial", design = bin_design)
      p11 <- predict(bin_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
      p00 <- predict(bin_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
      est.bart.sim.or[draw] <- (p11/(1-p11))/(p00/(1-p00))
      
      # kitchen sink approach
      
      Wmsm.out.ks.bart <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                        a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                                   data = dat1, method = "cbps",
                                   stabilize = TRUE)
      
      # estimate treatment effect
      
      bin_design <- svydesign(id=~1, weights = Wmsm.out.ks.bart$weights, data = dat2b)
      bin_mod <- svyglm(y ~ a_1*a_2,  family = "binomial", design = bin_design)
      p11 <- predict(bin_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
      p00 <- predict(bin_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
      est.bart.ks.or[draw] <- (p11/(1-p11))/(p00/(1-p00))
    }
  
    est.ps.sim.or <- unlist(est.ps.sim.or) # unlist to calculate mean and sd
    # results.it[1,14]<-mean(est.ps.sim.or)
    results.it[1,15]<-sd(est.ps.sim.or)
    results.it[1,16:17]<-c(results.it[1,14]-1.96*sd(est.ps.sim.or), results.it[1,14]+1.96*sd(est.ps.sim.or))
    
    est.ps.ks.or <- unlist(est.ps.ks.or)
    # results.it[1,18]<-mean(est.ps.ks.or)
    results.it[1,19]<-sd(est.ps.ks.or)
    results.it[1,20:21]<-c(results.it[1,18]-1.96*sd(est.ps.ks.or), results.it[1,18]+1.96*sd(est.ps.ks.or))
    
    est.cbps.sim.or <- unlist(est.cbps.sim.or)
    # results.it[1,22]<-mean(est.cbps.sim.or)
    results.it[1,23]<-sd(est.cbps.sim.or)
    results.it[1,24:25]<-c(results.it[1,22]-1.96*sd(est.cbps.sim.or), results.it[1,22]+1.96*sd(est.cbps.sim.or))
    
    est.cbps.ks.or <- unlist(est.cbps.ks.or)
    # results.it[1,26]<-mean(est.cbps.ks.or)
    results.it[1,27]<-sd(est.cbps.ks.or)
    results.it[1,28:29]<-c(results.it[1,26]-1.96*sd(est.cbps.ks.or), results.it[1,26]+1.96*sd(est.cbps.ks.or))
    
    est.bart.sim.or <- unlist(est.bart.sim.or)
    # results.it[1,30]<-mean(est.bart.sim.or)
    results.it[1,31]<-sd(est.bart.sim.or)
    results.it[1,32:33]<-c(results.it[1,30]-1.96*sd(est.bart.sim.or), results.it[1,30]+1.96*sd(est.bart.sim.or))
    
    est.bart.ks.or <- unlist(est.bart.ks.or)
    # results.it[1,34]<-mean(est.bart.ks.or)
    results.it[1,35]<-sd(est.bart.ks.or)
    results.it[1,36:37]<-c(results.it[1,34]-1.96*sd(est.bart.ks.or), results.it[1,34]+1.96*sd(est.bart.ks.or))
    
    # Baysian MSMs
    
    # JAGS
    
    # Bayesian inference 1. BMSM
    # first obtain MCMC sample for weights! from posterier distribution of treatment assignment parameters;
    jags.data<-list(w1=w1, w2=w2, a_1=a_1, a_1m= a_1, L1_1=L1_1,  L2_1=L2_1, 
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
    p1m<-matrix(NA, pdraws, ntot)
    p1c<-matrix(NA, pdraws, ntot)
    p2m<-matrix(NA, pdraws, ntot)
    p2c<-matrix(NA, pdraws, ntot)
    
    #calulating the MCMC weights;
    for (i2 in 1:(pdraws)){
      for (j2 in 1:(ntot)){
        
        p1m[i2,j2] <- expit(out.mcmc[j2,"bm10"])
        p1c[i2,j2]<- expit(out.mcmc[j2,"b10"] + out.mcmc[j2,"b11"]*w1[j2] + 
                             out.mcmc[j2,"b12"]*w1[j2] + out.mcmc[j2,"b13"]*L1_1[j2] + out.mcmc[j2,"b14"]*L2_1[j2])
        
        p2m[i2,j2] <- expit(out.mcmc[j2,"bm20"]+out.mcmc[j2,"bm21"]*a_1[j2])
        p2c[i2,j2] <- expit(out.mcmc[j2,"b20"]+out.mcmc[j2,"b21"]*w1[j2]+out.mcmc[j2,"b22"]*w2[j2]+
                              out.mcmc[j2,"b23"]*L1_1[j2]+out.mcmc[j2,"b24"]*L2_1[j2]+out.mcmc[j2,"b25"]*L1_2[j2]+
                              out.mcmc[j2,"b26"]*L2_2[j2]+out.mcmc[j2,"b27"]*a_1[j2])
        
        # # conditional treatment assignment model, v2;
        # a_2[i] ~ dbern(p2c[i])
        # logit(p2c[i]) <- b20 + b21*w1[i]+ b22*w2[i] + b23*L1_1[i] + b24*L2_1[i] + b25*L1_2[i] + b26*L2_2[i] + b27*a_1[i]
        # 
        # # marginal treatment assignment model, v2;
        # a_2m[i] ~ dbern(p2m[i])
        # logit(p2m[i]) <- bm20 + bm21*a_1[i]
        # 
        # # conditional treatment assignment model, v1;
        # a_1[i] ~ dbern(p1c[i])
        # logit(p1c[i]) <- b10 + b11*w1[i]+ b12*w2[i] + b13*L1_1[i] + b14*L2_1[i]
        # 
        # # marginal treatment assignment model, v1;
        # a_1m[i] ~ dbern(p1m[i])
        # logit(p1m[i]) <- bm10
        
        
        # exp_prob1[i2,j2] <- (exp(a_1[j2,1]*out.mcmc[i2,8]))/(1.0+exp(out.mcmc[i2,8]))
        # exp_prob2[i2,j2] <- exp_prob1[i2,j2]*(exp(z[j2,2]*(out.mcmc[i2,9]+out.mcmc[i2,10]*z[j2,1])))/(1.0+exp(out.mcmc[i2,9]+out.mcmc[i2,10]*z[j2,1]))
        
        # obs_prob1[i2,j2] <- (exp(z[j2,1]*(out.mcmc[i2,1] + out.mcmc[i2,2]*Lobs[j2,1]+out.mcmc[i2,3]*Lobs[j2,2])))/(1.0+exp(out.mcmc[i2,1] + out.mcmc[i2,2]*Lobs[j2,1]+out.mcmc[i2,3]*Lobs[j2,2]))
        # obs_prob2[i2,j2] <- obs_prob1[i2,j2]*(exp(z[j2,2]*(out.mcmc[i2,4]+out.mcmc[i2,5]*Lobs[j2,3]+out.mcmc[i2,6]*Lobs[j2,4]+out.mcmc[i2,7]*z[j2,1])))/(1.0+exp(out.mcmc[i2,4]+out.mcmc[i2,5]*Lobs[j2,3]+out.mcmc[i2,6]*Lobs[j2,4]+out.mcmc[i2,7]*z[j2,1]))
        
      }
      
      # if (i2 %% 50 == 0) {
      #   print(i2)
      # }
    }
    
    wmean <- colSums(p2m)/colSums(p2c)
    
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
      alpha <- as.numeric(rdirichlet(1, rep(1.0, length(dat2$y))))
      maxim <- optim(inits1,
                     fn=wloglik_binary,
                     Y=dat2$y,
                     A=cbind(1,dat2$a_1, dat2$a_2, dat2$a_1*dat2$a_2), #three mean parameters (intercept + coefficient for a_1 and coefficient for a_2);
                     weight=alpha*wmean,
                     control=list(fnscale=-1), method='BFGS', hessian=F)
      # need to estimate OR for logistic regression
      # (p11/(1-p11))/(p00/(1-p00))
      bootest2[j] <- (expit(maxim$par[1]+maxim$par[2]+maxim$par[3]+maxim$par[4])/(1-expit(maxim$par[1]+maxim$par[2]+maxim$par[3]+maxim$par[4])))/(expit(maxim$par[1])/(1-expit(maxim$par[1])))
      
      # OR we can report difference on the mean of Y between always treated and never treated;
      # bootest2[j] <- expit(maxim$par[1]+maxim$par[2]+maxim$par[3]+maxim$par[4]) - expit(maxim$par[1]) 
      
      # if (j %% 100 == 0) {
      #   print(j)
      # }
    }
    
    results.it[1,22] <- mean(bootest2)
    results.it[1,23] <- sd(bootest2)
    results.it[1,24:25] <- quantile(bootest2, probs=c(0.025,0.975))
    cbind(i,results.it)
  }
  
  # outfile <-"textcontsim"
  write.table(results_run, file = paste0(outfile,".txt"), row.names = FALSE,col.names = FALSE)
  
}

start_time <-Sys.time()
mysim.bin(paste0("Dec1","bin","run"), from=1, to=1000, ntot=1000, samplesize=10000, B=1000)
end_time <- Sys.time()
end_time - start_time

