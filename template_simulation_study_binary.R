###
#
###

library(tidyverse)
library(survey)
library(arm)
library(ltmle)
library(ipw) #MSM;
library(gfoRmula) #g-computation;
library(gtsummary)
library(SuperLearner)
library(WeightIt)
# library(confoundr)


#pre-define function;
expit <- function(x){
  x <- exp(x)/(exp(x)+1) 
  return(x)
}

mytrue <- function(a_1s = 1, a_2s = 1){
  #visit 1;
  w1sim <- rbinom(samplesize, 1, prob = 0.5) #50-50 male and female;
  w2sim <- rnorm(samplesize, mean = 12, sd = 4) #age;
  L1_1sim <- rbinom(samplesize, 1, prob = expit(-1 + 0.01*w1sim + 0.01*w2sim)) #a age related baseline binary clinical variable;
  L2_1sim <- rnorm(samplesize, mean = (0.01*w1sim + 0.01*w2sim), sd = 1) #a age & sex related baseline continuous clinical variable;
  
  #Visit 2, simulate all potential variables;
  L1_2sim <- rbinom(samplesize, 1, expit(-1 + 0.01*w1sim + 0.01*w2sim + 0.2*L1_1sim - 0.2*a_1s))
  L2_2sim <- rnorm(samplesize, mean = (L2_1sim - 0.01*w1sim + 0.01*w2sim - 0.2*a_1s), sd = 1)
  
  #Visit 3, simulate potential outcomes;
  # calc odds of y=1
  true_prob = expit(- 1 + log(0.9)*a_1s + log(0.85)*a_2s + log(1.01)*a_1s*a_2s 
                    + log(0.95)*w1sim + log(1.05)*w2sim + log(1.15)*L1_2sim + log(1.1)*L2_2sim)
  return(mean(true_prob))
}


# Simulation setup;
set.seed(123)
ntot=1000;
samplesize=10000;
B = 1000;
from = 1;
to = 1;

# placeholder matrix for outcome
results <- matrix(NA, 1000, 1+4*5?)

# true value;
# each method and each setting, est, se(est), low 95%CI, upper 95%CI;


  
for (i in from:to){ 

# data generation;
set.seed(i+123)

# Visit 1;
w1 <- rbinom(ntot, 1, prob = 0.5) #50-50 male and female (reference);
w2 <- rnorm(ntot, mean = 12, sd = 4) #age;
L1_1 <- rbinom(ntot, 1, prob = expit(-1 + 0.01*w1 + 0.01*w2)) #a age & sex related baseline binary clinical variable;
L2_1 <- rnorm(ntot, mean = (0.01*w1 + 0.01*w2), sd = 1) #a age & sex related baseline continuous clinical variable;

a_1 <- rbinom(ntot, 1, expit(-1 - 0.01*w1 + 0.01*w2 - 0.1*L1_1 + 0.1*L2_1)); #exposure;

# observational data;
L1_2 <- rbinom(ntot, 1, expit(-1 + 0.01*w1 + 0.01*w2 + 0.2*L1_1 - 0.2*a_1))
L2_2 <- rnorm(ntot, mean = (L2_1 - 0.01*w1 + 0.01*w2 - 0.2*a_1), sd = 1)
a_2 <- rbinom(ntot, 1, expit(-1 - 0.01*w1 + 0.01*w2 - 0.1*L1_2 + 0.1*L2_2 + a_1))
# table(a_1, a_2)

# end-of-study outcome;
y <- rbinom(ntot, 1, prob = expit(- 1 + log(0.9)*a_1 + log(0.85)*a_2 + log(1.01)*a_1*a_2 
                                  + log(0.95)*w1 + log(1.05)*w2 + log(1.15)*L1_2 + log(1.1)*L2_2))

# saving final data;
dat2 <- data.frame(w1, w2, L1_1, L2_1, a_1, L1_2, L2_2, a_2, y)


# convert prob into odds
# N0TE: plogis and expit are very close but not exactly the same
prob_11<-mytrue(a_1s = 1, a_2s = 1) #potential outcome Y(a_1=1, a_2=1);
prob_00<-mytrue(a_1s = 0, a_2s = 0)
# prob_10<-mytrue(a_1s = 1, a_2s = 0)
# prob_01<-mytrue(a_1s = 0, a_2s = 1)

# True absolute risk difference, -0.06414264; -0.06625605
# mean(prob_11) - mean(prob_00)
# True risk ratio (relative risk), 0.8425255; 0.8383083
# mean(prob_11)/mean(prob_00)
# True odds ratio, 0.7602477; 0.7537021

results[i,1]<-mean(prob_11/(1-prob_11))/mean(prob_00/(1-prob_00))



# gfoRmula package;
# no OR, need to do bootstrap
# first use the package it didnt automatically return OR 
# so we need to loop over the package to obtain the bt ci
# reached out to the developer -> update

# parametric gformula
est.gcomp.or <- rep(NA,B)

for (draw in 1:B){
  set.seed(draw)
  dat2b <- dat2[sample(1:ntot, size = ntot, replace = T),]
  dat2b <- tibble(dat2b)
  
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
  ymodel <- y ~ lag1_a + a + lag1_a*a + w1 + w2 + L1 + L2 + lag1_L1 + lag1_L2
  
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
  
  est.prob_00<-summary(gform_bin_eof)$result[1,4]
  est.prob_11<-summary(gform_bin_eof)$result[3,4]
  est.gcomp.or[draw] <- (est.prob_11/(1-est.prob_11))/(est.prob_00/(1-est.prob_00))
  
}

est.gcomp.or <- unlist(est.gcomp.or) # unlist to calculate mean and sd

results[i,2] <- mean(est.gcomp.or); # within the CI OR+1.96*OR_SE;
results[i,3] <- sd(est.gcomp.or)
results[i,4:5] <- quantile(est.gcomp.or, probs=c(0.025,0.975))

# ltmle package;
# has the option to specialy q and g
# OR we can use the SL (robust in terms of misspecification, non-linearity)
# can specify Q (cov and outcome) and G (treatment assignment model)
# more flexibility, more robust than the other two
# DO NOT work well in simulation study when sample size is small

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

results[i,6]<- out_tmle$effect.measures$OR$estimate
results[i,7]<- out_tmle$effect.measures$OR$std.dev
results[i,8:9] <- out_tmle$effect.measures$OR$CI


# correct Q form and correct G form
# a little bit underestimating
# doubly robust so large variance, wider CI

# tmle with superlearner + kitchen sink gform + Qform;
# great with large dataset
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
                       SL.library = "default", #with superlearner for continuous, SL.glm for binary, SL.random_forest, SL.gbm 
                       abar = list(c(1,1), c(0,0)),
                       estimate.time = FALSE)

out_tmle_s<-summary(tmle_model_SL, estimator="tmle")

results[i,10]<- out_tmle_s$effect.measures$OR$estimate
results[i,11]<- out_tmle_s$effect.measures$OR$std.dev
results[i,12:13] <- out_tmle_s$effect.measures$OR$CI

# ipw package;
# create long dataset
est.msm.or <- rep(NA, B)

for (draw in 1:B){
  #to be updated;
dat2 <- dat2 %>% 
  mutate(id = rep(1:1000),
         cum_a = a_1 + a_2)

dat2_new <- dat2 %>%
  pivot_longer(cols = -c(w1,w2,y,id,cum_a), 
               names_to = c("variable","visit"), 
               names_sep = "_", 
               values_to = "value") %>% 
  pivot_wider(names_from = variable, values_from = value) %>% 
  mutate(time = case_when(visit == 1 ~ 0,
                          visit == 2 ~ 1))

dat2_new$y[dat2_new$visit == 1] <- NA

dat2_new$visit <- as.numeric(dat2_new$visit)
dat2_new$w1 <- as.numeric(dat2_new$w1)

dat2_new$cum_a[dat2_new$time == 0] <- dat2$a_1
dat2_new$cum_a[dat2_new$time == 1] <- dat2$a_1+dat2$a_2

# calculate using ipwpoint

weights_v1 <- ipwpoint(
  exposure = a_1,
  family = "binomial",
  link = "logit",
  numerator = ~ 1,
  denominator = ~ w1 + w2 + L1_1 + L2_1,
  data = as.data.frame(dat2))

weights_v2 <- ipwpoint(
  exposure = a_2,
  family = "binomial",
  link = "logit",
  numerator = ~ a_1,
  denominator = ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1,
  data = as.data.frame(dat2))

w <- weights_v1$ipw.weights * weights_v2$ipw.weights

bin_design <- svydesign(id=~1, weights = ~ w, data = dat2)
bin_mod <- svyglm(y ~ as.factor(cum_a), family = "binomial", design = bin_design)
est.msm.or[draw]<-exp(coef(bin_mod))
}

est.msm.or <- unlist(est.msm.or) # unlist to calculate mean and sd
results[i,14]<-mean(est.msm.or)
results[i,15]<-sd(est.msm.or)
results[i,16:17]<-quantile(est.msm.or, probs=c(0.025,0.975))


# #calculating weights manually
# # IPT weights:
# 
# tmodel1 <- glm(a_1 ~ 1, family=binomial(link=logit), data = dat2)
# tmodel2 <- glm(a_2 ~ a_1, family=binomial(link=logit), data = dat2)
# tmodel2_package <- glm(a_2 ~ 1, family=binomial(link=logit), data = dat2)
# 
# smodel1 <- glm(a_1 ~ w1 + w2 + L1_1 + L2_1, 
#                family=binomial(link=logit), data = dat2)
# smodel2 <- glm(a_2 ~  w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1, 
#                family=binomial(link=logit), data = dat2) #kitchen sink approach
# smodel2_incorrect <- glm(a_2 ~  w1 + w2 +  L1_2 + L2_2,
#                          family=binomial(link=logit), data = dat2)
# 
# num <- predict(tmodel1, type = "response", data = dat2)* 
#   predict(tmodel2, type = "response", data = dat2)
# deno <- predict(smodel1, type = "response", data = dat2)* 
#   predict(smodel2, type = "response", data = dat2)
# deno_incorrect <- predict(smodel1, type = "response", data = dat2)*
#   predict(smodel2_incorrect, type = "response", data = dat2)
# 
# num_package <-predict(tmodel1, type = "response", data = dat2)* 
#   predict(tmodel2_package, type = "response", data = dat2)
# 
# deno_package<- predict(smodel1, type = "response", data = dat2)* 
#   predict(smodel2_incorrect, type = "response", data = dat2)
# 
# weights_package <- num_package/deno_package
# 
# correct_weights <- num/deno
# incorrect_weights <- num/deno_incorrect
# 
# library(survey)
# ds <- svydesign(id=~1, weights = ~ correct_weights, data = dat2)
# correct_fit <- svyglm(y ~ as.factor(cum_a), design = ds, family = "binomial")
# summary(correct_fit)
# exp(cbind(coef(correct_fit), confint(correct_fit)))
# 
# # calculate using the package's process
# 
# s <- split(w$ipw.weights, 1:2)
# v1 <- as.vector(s[[1]])
# v2 <- as.vector(s[[2]])
# m <- cbind(v1, v2)
# 
# ds_package <- svydesign(id=~1, weights = ~ m[,2], data = dat2)
# package_fit <- svyglm(y ~ as.factor(cum_a), design=ds_package, family = "binomial")
# summary(package_fit)
# exp(cbind(coef(package_fit), confint(package_fit)))
# 
# comparison <- tibble(cbind(m,  weights_package, correct_weights))
# View(comparison)
# 
# # Warning message:
# #   In eval(family$initialize) : non-integer #successes in a binomial glm!
# 
# #doesnt allow separate visit specification, require a long dataset

# WeightIt package

# bootstrap again here!;

# examine the initial imbalance at each time point and overall:
bal.tab(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
             a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
        data = dat2, stats = c("m", "ks"), thresholds = c(m = .05),
        which.time = .all)

# specify our weight models
Wmsm.out <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                             a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                        data = dat2, method = "ps",
                        stabilize = TRUE)

# take a look at the quality of the weights 
# summary(Wmsm.out)

# examine how well they perform with respect to covariate balance
# bal.tab(Wmsm.out, stats = c("m", "ks"), thresholds = c(m = .05),
        # which.time = .none)

# estimate treatment effect
# d.w.msm <- svydesign(~1, weights = Wmsm.out$weights,
                     # data = dat2)

cum.fit <- svyglm(y ~ as.factor(cum_a), design = d.w.msm, family = "binomial")
# Warning message: In eval(family$initialize) : non-integer #successes in a binomial glm!
summary(cum.fit)
exp(cbind(coef(cum.fit), confint(cum.fit)))
# exactly the same as manual ipw above

}



