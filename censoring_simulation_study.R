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

mysim.cens <- function(outfile, from=1, to=4, ntot=1000, samplesize=10000) {
  
  # Simulation setup;
  # ntot=1000;
  # samplesize=10000;
  # from = 1;
  # to = 3;
  
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
    truemean = (1 - 1*a_1s - 2*a_2s + 0.01*a_1s*a_2s - 0.001*w1sim + 0.01*w2sim + 0.1*L1_2sim + 0.1*L2_2sim)
    return(mean(truemean))
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
    # library(WeightIt)
    
    # data generation;
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
    
    # simulate c from a binary + logistic format;
    
    c <- rbinom(ntot, 1,
                prob = expit(-3 + 0.002*w1 + 0.001*w2 + 0.003*L1_1 + 0.005*L2_1 + 
                               0.01*a_1))
    
    dat4 <- data.frame(w1, w2, L1_1, L2_1, a_1, L1_2, L2_2, c, a_2, y)
    
    dat4$y[dat4$c == 1] <- NA
    dat4$a_2[dat4$c == 1] <- NA
    
    results.it[1,1] <- mytrue(a_1s = 1, a_2s = 1) - mytrue(a_1s = 0, a_2s = 0)
    
    # gfoRmula package;
    
    # create long table
    
    dat4_new <- dat4 %>%
      mutate(id = rep(1:1000)) %>% 
      pivot_longer(cols = -c(w1,w2,y,c,id), 
                   names_to = c("variable","visit"), 
                   names_sep = "_", 
                   values_to = "value") %>% 
      pivot_wider(names_from = variable, values_from = value) %>% 
      mutate(time = case_when(visit == 1 ~ 0,
                              visit == 2 ~ 1))
    
    dat4_new$y[dat4_new$visit == 1] <- NA
    
    id <- 'id'
    time_name <- 'time'
    covnames <- c("L1", "L2", "a")
    outcome_name <- 'y'
    covtypes <- c('binary', 'normal', 'binary')
    histories <- c(lagged)
    histvars <- list(c('a', 'L1', 'L2'))
    
    covparams <- list(covmodels = c(L1 ~ w1 + w2 + lag1_L1 + lag1_a,
                                    L2 ~ w1 + w2 + lag1_L2 + lag1_a,
                                    a ~ w1 + w2 + lag1_L1 + L1 + lag1_L2 + L2 + lag1_a))
    ymodel <- y ~ lag1_a*a + w1 + w2 + L1 + L2 + lag1_L1 + lag1_L2
    
    intvars <- list('a', 'a')
    interventions <- list(list(c(static, rep(0, 2))),
                          list(c(static, rep(1, 2))))
    int_descript <- c('Never treat', 'Always treat')
    
    censor_name <- 'c'
    censor_model <- c ~ w1 + w2 + lag1_L1 + lag1_L2 + lag1_a
    
    gform_censor <- gformula_continuous_eof(
      obs_data = dat4_new,
      id = id,
      time_name = time_name,
      covnames = covnames,
      outcome_name = outcome_name, 
      covtypes = c("binary", "normal", "binary"),
      covparams = covparams,  ymodel = ymodel,
      censor_name = censor_name, censor_model = censor_model,
      intvars = intvars, interventions = interventions,
      int_descript = int_descript, ref_int = 1,
      histories = c(lagged), histvars = list(c('a',"L1","L2")),
      boot_diag = TRUE, ## no difference in results
      basecovs = c("w1","w2"), 
      nsimul = 1000,
      nsamples = 1000, parallel = TRUE, ncores = 2,
      seed = 123)
    
    results.it[1, 2:5] <- unlist(gform_censor$result[3,12:15])
    
    # ltmle package;
    
    # move c column before y to prevent all columns after c being written as NA by the package
    # NA values are not permitted in data except after censoring or a survival event
    
    # change censoring indicator to satisfy package requirement
    dat4 <- dat4 %>% 
      mutate(c = case_when(c == 1 ~ "censored",
                           c == 0 ~ "uncensored"))
    
    dat4$c <- as.factor(dat4$c)
    
    # ltmle without superlearner + gform + Qform;
    tmle_model_noSL_c <- ltmle(dat4,
                             Anodes = c ("a_1","a_2"),
                             Cnodes = c("c"),
                             Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                             Ynodes = c("y"), 
                             survivalOutcome =FALSE,
                             Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a_1",
                                        y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a_1 + a_2 + a_1*a_2"),
                             gform = c("a_1 ~ w1 + w2 + L1_1 + L2_1",
                                       "c ~ w1 + w2 + L1_1 + L2_1 + a_1", # need to follow the order in data
                                       "a_2 ~ w1 + w2 + L1_1 + L2_1 + a_1 + L1_2 + L2_2"),
                             # gcomp = TRUE,
                             # iptw.only = FALSE,
                             # variance.method = "tmle",
                             # SL.library = "default", #with superlearner;
                             abar = list(c(1,1), c(0,0)),
                             estimate.time = FALSE)
    
    out_tmle <- summary(tmle_model_noSL_c, estimator="tmle")
    
    results.it[1,6]<- out_tmle$effect.measures$ATE$estimate
    results.it[1,7]<- out_tmle$effect.measures$ATE$std.dev
    results.it[1,8:9] <- out_tmle$effect.measures$ATE$CI
    
    # ltmle with superlearner + gform + Qform;
    tmle_model_SL_c <- ltmle(dat4,
                           Anodes = c ("a_1","a_2") ,
                           Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                           Ynodes = c("y"), 
                           Cnodes = c("c"),
                           survivalOutcome =FALSE,
                           Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a_1",
                                      y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a_1 + a_2 + a_1*a_2"),
                           gform = c("a_1 ~ w1 + w2 + L1_1 + L2_1",
                                     "c ~ w1 + w2 + L1_1 + L2_1 + a_1",
                                     "a_2 ~ w1 + w2 + L1_1 + L2_1 + a_1 + L1_2 + L2_2"),
                           # gcomp = TRUE,
                           # iptw.only = FALSE,
                           # variance.method = "tmle",
                           SL.library = "default", #with superlearner;
                           abar = list(c(1,1), c(0,0)),
                           estimate.time = FALSE)
    
    out_tmle_s<-summary(tmle_model_SL_c, estimator="tmle")
    
    results.it[1,10]<- out_tmle_s$effect.measures$ATE$estimate
    results.it[1,11]<- out_tmle_s$effect.measures$ATE$std.dev
    results.it[1,12:13] <- out_tmle_s$effect.measures$ATE$CI
    
    
    # MSM
    
    # simulation setup
    
    tmodel1 <- glm(a_1 ~ 1, family=binomial(link=logit), data = dat4)
    tmodel2 <- glm(a_2 ~ a_1, family=binomial(link=logit), data = dat4)
    
    smodel1 <- glm(a_1 ~ w1 + w2 + L1_1 + L2_1,
                   family=binomial(link=logit), data = dat4)
    smodel2.sim <- glm(a_2 ~  w1 + w2 + L1_2 + L2_2 + a_1,
                   family=binomial(link=logit), data = dat4) # simulation setup
    
    num <- predict(tmodel1, type = "response", newdata = dat4)*
      predict(tmodel2, type = "response", newdata = dat4)
    
    deno.sim <- predict(smodel1, type = "response", newdata = dat4)*
      predict(smodel2.sim, type = "response", newdata = dat4)
    
    weights.sim <- num/deno.sim
    
    # fit separate model to combine with weight
    cmodel <- glm(c ~ w1 + w2 + L1_1 + L2_1 + a_1, family=binomial(link=logit), data = dat4)
    csmodel <- glm(c ~ a_1, family=binomial(link=logit), data = dat4)
    
    cw <- predict(csmodel, type = "response", newdata = dat4)/predict(cmodel, type = "response", newdata = dat4)
    w.sim <- weights.sim * cw
    
    censor_sim_design <- svydesign(id=~1, weights = ~ w.sim, data = dat4)
    censor_sim_mod <- svyglm(y ~ a_1*a_2, design = censor_sim_design)
    p11 <- as.data.frame(predict(censor_sim_mod, newdata = data.frame(a_1=1, a_2=1)))[[1]]
    p11.se <- as.data.frame(predict(censor_sim_mod, newdata = data.frame(a_1=1, a_2=1)))[[2]]
    p00 <- as.data.frame(predict(censor_sim_mod, newdata = data.frame(a_1=0, a_2=0)))[[1]]
    p00.se <- as.data.frame(predict(censor_sim_mod, newdata = data.frame(a_1=0, a_2=0)))[[2]]
    results.it[1,14] <- p11 - p00
    results.it[1,15] <- sqrt(p11.se^2+p00.se^2)
    results.it[1,16:17] <- c(results.it[1,14]-1.96*sqrt(p11.se^2+p00.se^2),
                             results.it[1,14]+1.96*sqrt(p11.se^2+p00.se^2))
    
    
    # kitchen sink approach
    # 
    # tmodel1 <- glm(a_1 ~ 1, family=binomial(link=logit), data = dat4)
    # tmodel2 <- glm(a_2 ~ a_1, family=binomial(link=logit), data = dat4)
    # 
    # smodel1 <- glm(a_1 ~ w1 + w2 + L1_1 + L2_1,
    #                family=binomial(link=logit), data = dat4)
    smodel2 <- glm(a_2 ~  w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1,
                   family=binomial(link=logit), data = dat4) #kitchen sink approach

    # num <- predict(tmodel1, type = "response", newdata = dat4)*
    #   predict(tmodel2, type = "response", newdata = dat4)
    
    deno <- predict(smodel1, type = "response", newdata = dat4)*
      predict(smodel2, type = "response", newdata = dat4)
    
    weights <- num/deno
    
    # fit separate model to combine with weight
    # cmodel <- glm(c ~ w1 + w2 + L1_1 + L2_1 + a_1, family=binomial(link=logit), data = dat4)
    # csmodel <- glm(c ~ a_1, family=binomial(link=logit), data = dat4)
    # 
    # cw <- predict(csmodel, type = "response", newdata = dat4)/predict(cmodel, type = "response", newdata = dat4)
    w <- weights * cw
    
    censor_design <- svydesign(id=~1, weights = ~ w, data = dat4)
    censor_mod <- svyglm(y ~ a_1*a_2, design = censor_design)
    p11 <- as.data.frame(predict(censor_mod, newdata = data.frame(a_1=1, a_2=1)))[[1]]
    p11.se <- as.data.frame(predict(censor_mod, newdata = data.frame(a_1=1, a_2=1)))[[2]]
    p00 <- as.data.frame(predict(censor_mod, newdata = data.frame(a_1=0, a_2=0)))[[1]]
    p00.se <- as.data.frame(predict(censor_mod, newdata = data.frame(a_1=0, a_2=0)))[[2]]
    results.it[1,18] <- p11 - p00
    results.it[1,19] <- sqrt(p11.se^2+p00.se^2)
    results.it[1,20:21] <- c(results.it[1,18]-1.96*sqrt(p11.se^2+p00.se^2),
                             results.it[1,18]+1.96*sqrt(p11.se^2+p00.se^2))
    
    # 
    # est.sim <- rep(NA, B)
    # est.msm <- rep(NA, B)
    # 
    # 
    # for (draw in 1:B){
    #   set.seed(draw)
    #   dat4b <- dat4[sample(1:ntot, size = ntot, replace = T),]
    #   dat4b <- tibble(dat4b)
    #   
    #   # simulation setup
    #   
    #   tmodel1 <- glm(a_1 ~ 1, family=binomial(link=logit), data = dat4b)
    #   tmodel2 <- glm(a_2 ~ a_1, family=binomial(link=logit), data = dat4b)
    #   
    #   smodel1 <- glm(a_1 ~ w1 + w2 + L1_1 + L2_1,
    #                  family=binomial(link=logit), data = dat4b)
    #   smodel2.sim <- glm(a_2 ~  w1 + w2 + L1_2 + L2_2 + a_1,
    #                      family=binomial(link=logit), data = dat4b) # simulation setup
    #   
    #   num <- predict(tmodel1, type = "response", newdata = dat4b)*
    #     predict(tmodel2, type = "response", newdata = dat4b)
    #   
    #   deno.sim <- predict(smodel1, type = "response", newdata = dat4b)*
    #     predict(smodel2.sim, type = "response", newdata = dat4b)
    #   
    #   weights.sim <- num/deno.sim
    #   
    #   # fit separate model to combine with weight
    #   cmodel <- glm(c ~ w1 + w2 + L1_1 + L2_1 + a_1, family=binomial(link=logit), data = dat4b)
    #   csmodel <- glm(c ~ a_1, family=binomial(link=logit), data = dat4b)
    #   
    #   cw <- predict(csmodel, type = "response", newdata = dat4b)/predict(cmodel, type = "response", newdata = dat4b)
    #   w.sim <- weights.sim * cw
    #   
    #   censor_sim_design <- svydesign(id=~1, weights = ~ w.sim, data = dat4b)
    #   censor_sim_mod <- svyglm(y ~ a_1*a_2, design = censor_sim_design)
    #   p11 <- predict(censor_sim_mod, newdata = data.frame(a_1=1, a_2=1))[1]
    #   p00 <- predict(censor_sim_mod, newdata = data.frame(a_1=0, a_2=0))[1]
    #   est.sim[draw] <- p11-p00
    #   
    #   # kitchen sink approach
    # 
    #   smodel2 <- glm(a_2 ~  w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1,
    #                  family=binomial(link=logit), data = dat4b) #kitchen sink approach
    # 
    #   deno <- predict(smodel1, type = "response", newdata = dat4b)*
    #     predict(smodel2, type = "response", newdata = dat4b)
    #   
    #   weights <- num/deno
    #   
    #   w <- weights * cw
    #   
    #   censor_design <- svydesign(id=~1, weights = ~ w, data = dat4b)
    #   censor_mod <- svyglm(y ~ a_1*a_2, design = censor_design)
    #   p11 <- predict(censor_mod, newdata = data.frame(a_1=1, a_2=1))[1]
    #   p00 <- predict(censor_mod, newdata = data.frame(a_1=0, a_2=0))[1]
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
  write.table(results_run, file = paste0(outfile,".txt"), row.names = FALSE,col.names = FALSE)
  
}

start_time <-Sys.time()
mysim.cens(paste0("Apr30","contcens","run"), from=1, to=1000, ntot=1000, samplesize=10000)
end_time <- Sys.time()
end_time - start_time

