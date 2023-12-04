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

mysim.bin.cens <- function(outfile, from=1, to=4, ntot=1000, samplesize=10000, B=1000) {
  
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
    
    a_1 <- rbinom(ntot, 1, expit(-1 - 0.001*w1 + 0.01*w2 - 0.1*L1_1 + 0.1*L2_1)); #exposure;
    
    # observational data;
    L1_2 <- rbinom(ntot, 1, expit(-1 + 0.001*w1 + 0.01*w2 + 0.2*L1_1 - 0.2*a_1))
    L2_2 <- rnorm(ntot, mean = (L2_1 - 0.001*w1 + 0.01*w2 - 0.2*a_1), sd = 1)
    a_2 <- rbinom(ntot, 1, expit(-1 - 0.1*w1 + 0.01*w2 - 0.1*L1_2 + 0.01*L2_2 + 0.5*a_1))
    # table(a_1, a_2)
    
    # end-of-study outcome;
    y <- rbinom(ntot, 1, prob = expit(- 1 + log(0.9)*a_1 + log(0.85)*a_2 + log(1.01)*a_1*a_2 
                                      + log(0.95)*w1 + log(1.05)*w2 + log(1.15)*L1_2 + log(1.1)*L2_2))
    
    # simulate c from a binary + logistic format;
    
    c <- rbinom(ntot, 1,
                prob = expit(-3 + 0.002*w1 + 0.001*w2 + 0.003*L1_1 + 0.005*L2_1 + 
                               0.01*a_1))
    
    # saving final data;
    dat3 <- data.frame(w1, w2, L1_1, L2_1, a_1, L1_2, L2_2, c, a_2, y)
    
    dat3$y[dat3$c == 1] <- NA
    dat3$a_2[dat3$c == 1] <- NA
    
    # getting true OR from glm model using psuedo-population;
    ylong <- c(mytrue(a_1s = 0, a_2s = 0),mytrue(a_1s = 0, a_2s = 1),mytrue(a_1s = 1, a_2s = 0),mytrue(a_1s = 1, a_2s = 1))
    z1long <- c(rep(0, 2*samplesize), rep(1,2*samplesize))
    z2long <- c(rep(0, samplesize), rep(1,samplesize), rep(0, samplesize), rep(1, samplesize))
    
    truemodel <- glm(ylong ~ z1long*z2long)
    # summary(truemodel)
    prob_11 <-predict(truemodel, newdata = data.frame(z1long=1, z2long=1), type="response")
    prob_00 <-predict(truemodel, newdata = data.frame(z1long=0, z2long=0), type="response")
    
    results.it[1,1]<-(prob_11/(1-prob_11))/(prob_00/(1-prob_00))
    
    #gformula
    
    dat3_new <- dat3 %>%
      mutate(id = 1:1000) %>% 
      pivot_longer(cols = -c(w1,w2,y,c,id), 
                   names_to = c("variable","visit"), 
                   names_sep = "_", 
                   values_to = "value") %>% 
      pivot_wider(names_from = variable, values_from = value) %>% 
      mutate(time = case_when(visit == 1 ~ 0,
                              visit == 2 ~ 1)) # time has to start with 0
    
    dat3_new$y[dat3_new$visit == 1] <- NA
    
    id <- 'id' # must be named as id
    time_name <- 'time'
    covnames <- c("L1", "L2", "a")
    outcome_name <- 'y'
    covtypes <- c('binary', 'normal', 'binary')
    histories <- c(lagged)
    histvars <- list(c('a', 'L1', 'L2'))
    
    # have to specify each model for time dep cov, trt and outcome
    # need to parametrically specify all the time dep component
    covparams <- list(covmodels = c(L1 ~ w1 + w2 + lag1_L1 + lag1_a,
                                    L2 ~ w1 + w2 + lag1_L2 + lag1_a,
                                    a ~ w1 + w2 + lag1_L1 + L1 + lag1_L2 + L2 + lag1_a))
    ymodel <- y ~  lag1_a*a + w1 + w2 + L1 + L2 + lag1_L1 + lag1_L2
    
    censor_name <- 'c'
    censor_model <- c ~ w1 + w2 + lag1_L1 + lag1_L2 + lag1_a
    
    intvars <- list('a', 'a')
    
    interventions <- list(list(c(static, rep(0, 2))),
                          list(c(static, rep(1, 2))))
    int_descript <- c('Never treat', 'Always treat')
    
    gform_bin_eof <- gformula_binary_eof(
      obs_data = dat3_new,
      id = id,
      time_name = time_name,
      covnames = covnames,
      outcome_name = outcome_name, 
      covtypes = c("binary", "normal", "binary"),
      covparams = covparams,  
      ymodel = ymodel,
      censor_name = censor_name, censor_model = censor_model,
      intvars = intvars, interventions = interventions,
      int_descript = int_descript, ref_int = 1,
      histories = c(lagged), histvars = list(c('a',"L1","L2")),
      basecovs = c("w1","w2"),
      seed=123)
    
    est.prob_00<- unlist(summary(gform_bin_eof)$result[2,4])
    est.prob_11<- unlist(summary(gform_bin_eof)$result[3,4])
    results.it[1,2] <- (est.prob_11/(1-est.prob_11))/(est.prob_00/(1-est.prob_00))
    
    # gfoRmula package;
    est.gcomp.or <- rep(NA,B)
    
    for (draw in 1:B){
      set.seed(draw)
      dat3b <- dat3[sample(1:ntot, size = ntot, replace = T),]
      dat3b <- tibble(dat3b)
      
      # creating long format data;
      dat3_new_b <- dat3b %>%
        mutate(id = 1:1000) %>% 
        pivot_longer(cols = -c(w1,w2,y,c,id), 
                     names_to = c("variable","visit"), 
                     names_sep = "_", 
                     values_to = "value") %>% 
        pivot_wider(names_from = variable, values_from = value) %>% 
        mutate(time = case_when(visit == 1 ~ 0,
                                visit == 2 ~ 1)) # time has to start with 0
      
      dat3_new_b$y[dat3_new_b$visit == 1] <- NA
      
      id <- 'id' # must be named as id
      time_name <- 'time'
      covnames <- c("L1", "L2", "a")
      outcome_name <- 'y'
      covtypes <- c('binary', 'normal', 'binary')
      histories <- c(lagged)
      histvars <- list(c('a', 'L1', 'L2'))
      
      # have to specify each model for time dep cov, trt and outcome
      # need to parametrically specify all the time dep component
      covparams <- list(covmodels = c(L1 ~ w1 + w2 + lag1_L1 + lag1_a,
                                      L2 ~ w1 + w2 + lag1_L2 + lag1_a,
                                      a ~ w1 + w2 + lag1_L1 + L1 + lag1_L2 + L2 + lag1_a))
      ymodel <- y ~ lag1_a + a + lag1_a*a + w1 + w2 + L1 + L2 + lag1_L1 + lag1_L2
      
      censor_name <- 'c'
      censor_model <- c ~ w1 + w2 + lag1_L1 + lag1_L2 + lag1_a
      
      intvars <- list('a', 'a')
      
      interventions <- list(list(c(static, rep(0, 2))),
                            list(c(static, rep(1, 2))))
      int_descript <- c('Never treat', 'Always treat')
      
      gform_bin_eof <- gformula_binary_eof(
        obs_data = dat3_new_b,
        id = id,
        time_name = time_name,
        covnames = covnames,
        outcome_name = outcome_name, 
        covtypes = c("binary", "normal", "binary"),
        covparams = covparams,  
        ymodel = ymodel,
        censor_name = censor_name, censor_model = censor_model,
        intvars = intvars, interventions = interventions,
        int_descript = int_descript, ref_int = 1,
        histories = c(lagged), histvars = list(c('a',"L1","L2")),
        basecovs = c("w1","w2"),
        seed=123)
      
      est.prob_00<-summary(gform_bin_eof)$result[1,4]
      est.prob_11<-summary(gform_bin_eof)$result[3,4]
      est.gcomp.or[draw] <- (est.prob_11/(1-est.prob_11))/(est.prob_00/(1-est.prob_00))
      
    }
    
    est.gcomp.or <- unlist(est.gcomp.or) # unlist to calculate mean and sd
    
    # results.it[1,2] <- mean(est.gcomp.or); # within the CI OR+1.96*OR_SE;
    results.it[1,3] <- sd(est.gcomp.or)
    results.it[1,4:5] <- c(results.it[1,2]-1.96*sd(est.gcomp.or), results.it[1,2]+1.96*sd(est.gcomp.or))
    
    # ltmle package;
    
    dat3 <- dat3 %>% 
      mutate(c = case_when(c == 1 ~ "censored",
                           c == 0 ~ "uncensored"))
    
    dat3$c <- as.factor(dat3$c)
    
    tmle_model_noSL <- ltmle(dat3,
                             Anodes = c ("a_1","a_2") ,
                             Cnodes = c("c"),
                             Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                             Ynodes = c("y"), 
                             survivalOutcome =FALSE,
                             Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a_1",
                                        y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a_1 + a_2 + a_1*a_2"),
                             gform = c("a_1 ~ w1 + w2 + L1_1 + L2_1",
                                       "c ~ w1 + w2 + L1_1 + L2_1 + a_1",
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
    tmle_model_SL <- ltmle(dat3,
                           Anodes = c ("a_1","a_2") ,
                           Cnodes = c("c"),
                           Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                           Ynodes = c("y"), 
                           survivalOutcome =FALSE,
                           Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a_1",
                                      y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a_1 + a_2 + a_1*a_2"),
                           gform = c("a_1 ~ w1 + w2 + L1_1 + L2_1",
                                     "c ~ w1 + w2 + L1_1 + L2_1 + a_1",
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
    
    # MSM
    
    # simulation setup
    
    tmodel1 <- glm(a_1 ~ 1, family=binomial(link=logit), data = dat3)
    tmodel2 <- glm(a_2 ~ a_1, family=binomial(link=logit), data = dat3)
    
    smodel1 <- glm(a_1 ~ w1 + w2 + L1_1 + L2_1,
                   family=binomial(link=logit), data = dat3)
    smodel2.sim <- glm(a_2 ~  w1 + w2 + L1_2 + L2_2 + a_1,
                       family=binomial(link=logit), data = dat3) # simulation setup
    
    num <- predict(tmodel1, type = "response", newdata = dat3)*
      predict(tmodel2, type = "response", newdata = dat3)
    
    deno.sim <- predict(smodel1, type = "response", newdata = dat3)*
      predict(smodel2.sim, type = "response", newdata = dat3)
    
    weights.sim <- num/deno.sim
    
    # fit separate model to combine with weight
    cmodel <- glm(c ~ w1 + w2 + L1_1 + L2_1 + a_1, family=binomial(link=logit), data = dat3)
    csmodel <- glm(c ~ a_1, family=binomial(link=logit), data = dat3)
    cw <- predict(csmodel, type = "response", newdata = dat3)/predict(cmodel, type = "response", newdata = dat3)
    
    w.sim <- weights.sim * cw
    
    bin_sim_design <- svydesign(id=~1, weights = ~ w.sim, data = dat3)
    bin_sim_mod <- svyglm(y ~ a_1*a_2, family = "binomial", design = bin_sim_design)
    p11 <- predict(bin_sim_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
    p00 <- predict(bin_sim_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
    results.it[1,14] <- (p11/(1-p11))/(p00/(1-p00))
    
    # kitchen sink approach
    
    # tmodel1 <- glm(a_1 ~ 1, family=binomial(link=logit), data = dat3)
    # tmodel2 <- glm(a_2 ~ a_1, family=binomial(link=logit), data = dat3)
    # 
    # smodel1 <- glm(a_1 ~ w1 + w2 + L1_1 + L2_1,
    #                family=binomial(link=logit), data = dat3)
    smodel2 <- glm(a_2 ~  w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1,
                   family=binomial(link=logit), data = dat3) #kitchen sink approach
    # 
    # num <- predict(tmodel1, type = "response", newdata = dat3)*
    #   predict(tmodel2, type = "response", newdata = dat3)
    
    deno <- predict(smodel1, type = "response", newdata = dat3)*
      predict(smodel2, type = "response", newdata = dat3)
    
    weights <- num/deno
    
    # fit separate model to combine with weight
    # cmodel <- glm(c ~ w1 + w2 + L1_1 + L2_1 + a_1, family=binomial(link=logit), data = dat3)
    # csmodel <- glm(c ~ a_1, family=binomial(link=logit), data = dat3)
    # cw <- predict(csmodel, type = "response", newdata = dat3)/predict(cmodel, type = "response", newdata = dat3)
    
    w <- weights * cw
    
    bin_sim_design <- svydesign(id=~1, weights = ~ w, data = dat3)
    bin_sim_mod <- svyglm(y ~ a_1*a_2, family = "binomial", design = bin_sim_design)
    p11 <- predict(bin_sim_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
    p00 <- predict(bin_sim_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
    results.it[1,18] <- (p11/(1-p11))/(p00/(1-p00))
   
    
    est.sim.or <- rep(NA, B)
    est.msm.or <- rep(NA, B)
    

    for (draw in 1:B){
      set.seed(draw)
      dat3b <- dat3[sample(1:ntot, size = ntot, replace = T),]
      dat3b <- tibble(dat3b)
      
      # simulation setup
      
      tmodel1 <- glm(a_1 ~ 1, family=binomial(link=logit), data = dat3b)
      tmodel2 <- glm(a_2 ~ a_1, family=binomial(link=logit), data = dat3b)
      
      smodel1 <- glm(a_1 ~ w1 + w2 + L1_1 + L2_1,
                     family=binomial(link=logit), data = dat3b)
      smodel2.sim <- glm(a_2 ~  w1 + w2 + L1_2 + L2_2 + a_1,
                         family=binomial(link=logit), data = dat3b) # simulation setup
      
      num <- predict(tmodel1, type = "response", newdata = dat3b)*
        predict(tmodel2, type = "response", newdata = dat3b)
      
      deno.sim <- predict(smodel1, type = "response", newdata = dat3b)*
        predict(smodel2.sim, type = "response", newdata = dat3b)
      
      weights.sim <- num/deno.sim
      
      # fit separate model to combine with weight
      cmodel <- glm(c ~ w1 + w2 + L1_1 + L2_1 + a_1, family=binomial(link=logit), data = dat3b)
      csmodel <- glm(c ~ a_1, family=binomial(link=logit), data = dat3b)
      
      cw <- predict(csmodel, type = "response", newdata = dat3b)/predict(cmodel, type = "response", newdata = dat3b)
      w.sim <- weights.sim * cw
      
      bin_sim_design <- svydesign(id=~1, weights = ~ w.sim, data = dat3b)
      bin_sim_mod <- svyglm(y ~ a_1*a_2, family = "binomial", design = bin_sim_design)
      p11 <- predict(bin_sim_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
      p00 <- predict(bin_sim_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
      est.sim.or[draw] <- (p11/(1-p11))/(p00/(1-p00))
      
      # kitchen sink approach
      
      tmodel1 <- glm(a_1 ~ 1, family=binomial(link=logit), data = dat3b)
      tmodel2 <- glm(a_2 ~ a_1, family=binomial(link=logit), data = dat3b)
      
      smodel1 <- glm(a_1 ~ w1 + w2 + L1_1 + L2_1,
                     family=binomial(link=logit), data = dat3b)
      smodel2 <- glm(a_2 ~  w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1,
                     family=binomial(link=logit), data = dat3b) #kitchen sink approach
      
      num <- predict(tmodel1, type = "response", newdata = dat3b)*
        predict(tmodel2, type = "response", newdata = dat3b)
      
      deno <- predict(smodel1, type = "response", newdata = dat3b)*
        predict(smodel2, type = "response", newdata = dat3b)
      
      weights <- num/deno
      
      # fit separate model to combine with weight
      cmodel <- glm(c ~ w1 + w2 + L1_1 + L2_1 + a_1, family=binomial(link=logit), data = dat3b)
      csmodel <- glm(c ~ a_1, family=binomial(link=logit), data = dat3b)
      
      cw <- predict(csmodel, type = "response", newdata = dat3b)/predict(cmodel, type = "response", newdata = dat3b)
      w <- weights * cw
      
      bin_design <- svydesign(id=~1, weights = ~ w, data = dat3b)
      bin_mod <- svyglm(y ~ a_1*a_2, family = "binomial", design = bin_design)
      p11 <- predict(bin_mod, newdata = data.frame(a_1=1, a_2=1), type = "response")[1]
      p00 <- predict(bin_mod, newdata = data.frame(a_1=0, a_2=0), type = "response")[1]
      est.msm.or[draw] <- (p11/(1-p11))/(p00/(1-p00))
      
    }
    
    est.sim.or <- unlist(est.sim.or)
    # results.it[1,14]<-mean(est.sim.or)
    results.it[1,15]<-sd(est.sim.or)
    results.it[1,16:17]<-c(results.it[1,14]-1.96*sd(est.sim.or), results.it[1,14]+1.96*sd(est.sim.or))
    
    est.msm.or <- unlist(est.msm.or) # unlist to calculate mean and sd
    # results.it[1,18]<-mean(est.msm.or)
    results.it[1,19]<-sd(est.msm.or)
    results.it[1,20:21]<-c(results.it[1,18]-1.96*sd(est.msm.or), results.it[1,18]+1.96*sd(est.msm.or))
    
    cbind(i,results.it)
  }
  
  # outfile <-"textcontsim"
  write.table(results_run, file = paste0(outfile,".txt"), row.names = FALSE,col.names = FALSE)
  
}

start_time <-Sys.time()
mysim.bin.cens(paste0("Dec2","bin","cens","run"), from=1, to=1000, ntot=1000, samplesize=10000, B=1000)
end_time <- Sys.time()
end_time - start_time


