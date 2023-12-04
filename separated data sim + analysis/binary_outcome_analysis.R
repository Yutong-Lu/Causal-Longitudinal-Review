##
#
##

library(tidyverse)
library(survey)
library(arm)
library(ltmle)
library(ipw) #MSM;
library(gfoRmula) #g-computation;
library(gtsummary)
library(SuperLearner)
library(WeightIt)
library(confoundr)
options(scipen=999, warn = -1)

dat2 <- read_csv("datasets/binary_outcome_data.csv")


# gfoRmula package;
# no OR, need to do bootstrap
# first use the package it didnt automatically return OR 
# so we need to loop over the package to obtain the bt ci
# reached out to the developer -> update

# parametric gformula

n = length(dat2$y)
B = 1000;
est.or <- rep(NA,B)

for (draw in 1:B){
  set.seed(draw)
  
  dat2b <- dat2[sample(1:n, size = n, replace = T),]
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
  est.or[draw] <- (est.prob_11/(1-est.prob_11))/(est.prob_00/(1-est.prob_00))

}

est.or <- unlist(est.or) # unlist to calculate mean and sd

OR <- mean(est.or); # within the CI OR+1.96*OR_SE
OR_SE <- sd(est.or)

summary(gform_bin_eof)
OR
OR_SE
# 95% CI
OR+1.96*OR_SE
OR-1.96*OR_SE

# or;
# (0.3019694/(1-0.3019694))/(0.3884557/(1-0.3884557)); 0.6810436;

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

summary(tmle_model_noSL, estimator="tmle")
# correct Q form and correct G form
# a little bit underestimating
# doubly robust so large variance, wider CI

# tmle with superlearner + kitchen sink gform + Qform;
# great with large dataset
set.seed(123)
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

summary(tmle_model_SL, estimator="tmle")

#tmle package for gcomputation + kitchen sink gform + correct Qform;
tmle_model_gcomp <- ltmle(dat2,
                                   Anodes = c ("a_1","a_2") ,
                                   Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                                   Ynodes = c("y"), 
                                   survivalOutcome =FALSE,
                                   Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a_1",
                                              y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a_1 + a_2 + a_1*a_2"),
                                   gform = c("a_1 ~ w1 + w2 + L1_1 + L2_1",
                                             "a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1"),
                                   gcomp = TRUE,
                                   # iptw.only = FALSE,
                                   # variance.method = "tmle",
                                   abar = list(c(1,1), c(0,0)),
                                   estimate.time = FALSE)

summary(tmle_model_gcomp)
# tmle package can also run gcomp but recommend using gformula 
# (people should use the dedicated package)

# ipw package;

# create long dataset
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

w <- ipwtm(
  exposure = a,
  family = "binomial",
  link = "logit",
  numerator = ~ w1 + w2,
  denominator = ~ w1 + w2 + L1 + L2,
  id = id,
  timevar = visit,
  type = "first",
  data = as.data.frame(dat2_new))

summary(w)

#plot inverse probability weights
graphics.off()
ipwplot(weights = w$ipw.weights, timevar = dat2_new$visit,
        binwidth = 1, ylim = c(-1.5, 1.5), main = "Stabilized inverse probability weights")

# need a design
ds_weights_directly_from_package <- svydesign(id=~1, weights = ~ w$ipw.weights, data = dat2_new)
summary(svyglm(y ~ as.factor(cum_a), data = dat2_new, family = "binomial", 
               weights = w$ipw.weights, design=ds_weights_directly_from_package))
# log OR to OR
exp(coef(summary(svyglm(y ~ as.factor(cum_a), data = dat2_new, family = "binomial", 
                        weights = w$ipw.weights, design=ds_weights_directly_from_package))))

# model using ipwpoint
# weight with time dep treatment model -> calculate by visit and combine
# break up the weight calculation using the package
# no obvious adv using the package vs manual calculation
# ppl need to know the method a little bit to do it manually

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
summary(bin_mod)
confint(bin_mod)
exp(cbind(coef(bin_mod), confint(bin_mod)))

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
summary(Wmsm.out)

# examine how well they perform with respect to covariate balance
bal.tab(Wmsm.out, stats = c("m", "ks"), thresholds = c(m = .05),
        which.time = .none)

# estimate treatment effect
d.w.msm <- svydesign(~1, weights = Wmsm.out$weights,
                     data = dat2)

cum.fit <- svyglm(y ~ as.factor(cum_a), design = d.w.msm, family = "binomial")
# Warning message: In eval(family$initialize) : non-integer #successes in a binomial glm!
summary(cum.fit)
exp(cbind(coef(cum.fit), confint(cum.fit)))
# exactly the same as manual ipw above

# Confoundr package

# DEMONSTRATION OF DIAGNOSTIC 1

# create the appropriate history strata
dat2.history <- makehistory.one(
  input=dat2,
  id="id",
  exposure="a",
  times=c(1,2),
  name.history="h"
)

# create a tidy dataframe (long data)
dat2.tidy <- confoundr::lengthen(
  input=dat2.history,
  id="id",
  diagnostic=1,
  censoring="no",
  exposure="a",
  temporal.covariate=c("L1", "L2"),
  static.covariate = c("w1", "w2"),
  times.exposure=c(1,2), 
  times.covariate=c(1,2),
  history="h"
)

# create a covariate balance table
baltbl2 <- confoundr::balance(
  input=dat2.tidy,
  diagnostic=1,
  approach="none",
  censoring="no",
  scope="recent",
  recency=0,
  exposure="a",
  history="h",
  times.exposure=c(1,2),
  times.covariate=c(1,2),
  sd.ref="yes"
)

# create a trellised covariate balance plot
balplot2 <- confoundr::makeplot(
  input=baltbl2,
  diagnostic=1,
  approach="none",
  censoring="no",
  scope="recent",
  label.exposure="Treatment",
  label.covariate="Covariate"
)

balplot2

# DEMONSTRATION OF DIAGNOSTIC 3

dat2.history$w_1 <- weights_v1$ipw.weights
dat2.history$w_2 <- weights_v2$ipw.weights

dat2.history3 <- confoundr::lengthen(
  input=dat2.history,
  id="id",
  diagnostic=3,
  censoring="no",
  exposure="a",
  temporal.covariate=c("L1", "L2"),
  static.covariate = c("w1", "w2"),
  times.exposure=c(1,2), 
  times.covariate=c(1,2),
  history="h",
  weight.exposure="w"
)

baltbl2_diag3 <- confoundr::balance(
  input=dat2.history3,
  diagnostic=3,
  approach="weight",
  censoring="no",
  scope="all",
  exposure="a",
  history="h",
  times.exposure=c(1,2),
  times.covariate=c(1,2),
  sd.ref="yes",
  weight.exposure="w"
)

balplot2_diag3 <- confoundr::makeplot(
  input=baltbl2_diag3,
  diagnostic=3,
  approach="weight",
  censoring="no",
  scope="all",
  label.exposure="Study Dropout",
  label.covariate="Covariate"
)

balplot2_diag3

# # End
