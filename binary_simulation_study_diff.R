###
#
###
setwd("/Users/lorrainelu/Documents/GitHub/longitudinal_causal_analysis")
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
  
    
    
    # data generation;
    set.seed(i+123)
    
    results.it <- matrix(NA, 1, 21)
  
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
    
    results.it[1,1]<-(prob_11-prob_00)
    
      
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
    results.it[1,2] <- (est.prob_11-est.prob_00)
  
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
      est.gcomp.or[draw] <- (est.prob_11-est.prob_00)
  
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
  
    results.it[1,6]<- out_tmle$effect.measures$ATE$estimate
    results.it[1,7]<- out_tmle$effect.measures$ATE$std.dev
    results.it[1,8:9] <- out_tmle$effect.measures$ATE$CI
  
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
  
    results.it[1,10]<- out_tmle_s$effect.measures$ATE$estimate
    results.it[1,11]<- out_tmle_s$effect.measures$ATE$std.dev
    results.it[1,12:13] <- out_tmle_s$effect.measures$ATE$CI
  
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
    results.it[1,14] <- (p11-p00)
  
    # kitchen sink
  
    Wmsm.out.sim <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                     a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                                data = dat2, method = "ps",
                                stabilize = TRUE)
  
    bin_design <- svydesign(id=~1, weights = Wmsm.out.sim$weights, data = dat2)
    bin_mod <- svyglm(y ~ a_1*a_2, family = "binomial", design = bin_design)
    p11 <- predict(bin_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
    p00 <- predict(bin_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
    results.it[1,18] <- (p11-p00)
  
    est.msm.or <- rep(NA, B)
    est.weightit.or <- rep(NA, B)
  
    for (draw in 1:B){
      set.seed(draw)
      dat2b <- dat2[sample(1:ntot, size = ntot, replace = T),]
      dat2b <- tibble(dat2b)
  
      # dat2b <- dat2b %>%
      #   mutate(id = rep(1:1000),
      #          cum_a = a_1 + a_2)
  
      # calculate using simulation setup
    
      Wmsm.out.sim <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                       a_2 ~ w1 + w2 + L1_2 + L2_2 + a_1),
                                  data = dat2b, method = "ps",
                                  stabilize = TRUE)
    
      bin_design <- svydesign(id=~1, weights = Wmsm.out.sim$weights, data = dat2b)
      bin_mod <- svyglm(y ~ a_1*a_2, family = "binomial", design = bin_design)
      p11 <- predict(bin_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
      p00 <- predict(bin_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
      est.msm.or[draw] <- (p11-p00)
    
      # calculate using kitchen sink approach
    
      Wmsm.out <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                   a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                              data = dat2b, method = "ps",
                              stabilize = TRUE)
    
      d.w.msm <- svydesign(~1, weights = Wmsm.out$weights, data = dat2b)
      cum.fit <- svyglm(y ~ a_1*a_2, design = d.w.msm, family = "binomial")
      p11 <- predict(cum.fit, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
      p00 <- predict(cum.fit, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
      est.weightit.or[draw] <- (p11-p00)
    }
  
    est.msm.or <- unlist(est.msm.or) # unlist to calculate mean and sd
    # results.it[1,14]<-mean(est.msm.or)
    results.it[1,15]<-sd(est.msm.or)
    results.it[1,16:17]<-c(results.it[1,14]-1.96*sd(est.msm.or), results.it[1,14]+1.96*sd(est.msm.or))
    
    est.weightit.or <- unlist(est.weightit.or)
    # results.it[1,18]<-mean(est.weightit.or)
    results.it[1,19]<-sd(est.weightit.or)
    results.it[1,20:21]<-c(results.it[1,18]-1.96*sd(est.weightit.or), results.it[1,18]+1.96*sd(est.weightit.or))
  
    cbind(i,results.it)
  }
  
  # outfile <-"textcontsim"
  write.table(results_run, file = paste0(outfile,".txt"), row.names = FALSE,col.names = FALSE)
  
}

start_time <-Sys.time()
mysim.bin(paste0("Jan30","bin","run"), from=1, to=100, ntot=1000, samplesize=10000, B=1000)
end_time <- Sys.time()
end_time - start_time

