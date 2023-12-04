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
options(scipen=999)

dat4 <- read_csv("datasets/censoring_data.csv")

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
ymodel <- y ~ lag1_a + a + lag1_a*a + w1 + w2 + L1 + L2 + lag1_L1 + lag1_L2

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

summary(gform_censor)

# ltmle package;

# move c column before y to prevent all columns after c being written as NA by the package
# NA values are not permitted in data except after censoring or a survival event

# change censoring indicator to satisfy package requirement
dat4 <- dat4 %>% 
  mutate(c = case_when(c == 1 ~ "censored",
                       c == 0 ~ "uncensored"))

dat4$c <- as.factor(dat4$c)

# ltmle without superlearner + gform + Qform;
tmle_model_noSL <- ltmle(dat4,
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

summary(tmle_model_noSL, estimator="tmle")

# ltmle with superlearner + gform + Qform;
set.seed(123)
tmle_model_SL <- ltmle(dat4,
                    Anodes = c ("a_1","a_2") ,
                    Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                    Ynodes = c("y"), 
                    Cnodes = c("c"),
                    survivalOutcome =FALSE,
                    Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a_1",
                               y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a_1 + a_2 + a_1*a_2"),
                    gform = c("a_1 ~ w1 + w2 + L1_1 + L2_1",
                              "c ~ w1 + w2 + L1_1 + L2_1 + a_1", # need to follow the order in data
                              "a_2 ~ w1 + w2 + L1_1 + L2_1 + a_1 + L1_2 + L2_2"),
                    # gcomp = TRUE,
                    # iptw.only = FALSE,
                    # variance.method = "tmle",
                    SL.library = "default", #with superlearner;
                    abar = list(c(1,1), c(0,0)),
                    estimate.time = FALSE)

summary(tmle_model_SL, estimator="tmle")

#ltmle package for gcomputation + gform + Qform;
tmle_model_gcomp <- ltmle(dat4,
                    Anodes = c ("a_1","a_2") ,
                    Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                    Ynodes = c("y"), 
                    Cnodes = c("c"),
                    survivalOutcome =FALSE,
                    Qform = c(L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a_1",
                               y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a_1 + a_2 + a_1*a_2"),
                    gform = c("a_1 ~ w1 + w2 + L1_1 + L2_1",
                              "c ~ w1 + w2 + L1_1 + L2_1 + a_1", # need to follow the order in data
                              "a_2 ~ w1 + w2 + L1_1 + L2_1 + a_1 + L1_2 + L2_2"),
                    gcomp = TRUE,
                    # iptw.only = FALSE,
                    # variance.method = "tmle",
                    abar = list(c(1,1), c(0,0)),
                    estimate.time = FALSE)

summary(tmle_model_gcomp)

# iptw

dat4 <- dat4 %>% 
  mutate(id = rep(1:1000),
         cum_a = a_1 + a_2)

dat4_new <- dat4 %>%
  pivot_longer(cols = -c(w1,w2,y,id,cum_a,c), 
               names_to = c("variable","visit"), 
               names_sep = "_", 
               values_to = "value") %>% 
  pivot_wider(names_from = variable, values_from = value) %>% 
  mutate(time = case_when(visit == 1 ~ 0,
                          visit == 2 ~ 1))

dat4_new$y[dat4_new$visit == 1] <- NA

dat4_new$visit <- as.numeric(dat4_new$visit)
dat4_new$w1 <- as.numeric(dat4_new$w1)

dat4_new$cum_a[dat4_new$time == 0] <- dat4$a_1
dat4_new$cum_a[dat4_new$time == 1] <- dat4$a_1+dat4$a_2

# w <- ipwtm(
#   exposure = a,
#   family = "binomial",
#   link = "logit",
#   numerator = ~ 1,
#   denominator = ~ w1 + w2 + L1 + L2,
#   id = id,
#   timevar = visit,
#   type = "all",
#   data = as.data.frame(dat4_new))

# use wide data set and calculate at each visit

weights_v1 <- ipwpoint(
  exposure = a_1,
  family = "binomial",
  link = "logit",
  numerator = ~ 1,
  denominator = ~ w1 + w2 + L1_1 + L2_1,
  data = as.data.frame(dat4))

weights_v2 <- ipwpoint(
  exposure = a_2,
  family = "binomial",
  link = "logit",
  numerator = ~ a_1,
  denominator = ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1,
  data = as.data.frame(dat4)) # doesnt work with NAs in a_2

# fit separate model to combine with weight
cmodel <- glm(c ~ w1 + w2 + L1_1 + L2_1 + a_1, family=binomial(link=logit), data = dat4)
csmodel <- glm(c ~ a_1, family=binomial(link=logit), data = dat4)

cw <- predict(csmodel, type = "response", newdata = dat4)/predict(cmodel, type = "response", newdata = dat4)
w <- weights_v1$ipw.weights * weights_v2$ipw.weights * cw

censor_design <- svydesign(id=~1, weights = ~ w, data = dat4)
censor_mod <- svyglm(y ~ as.factor(cum_a), design = censor_design)
summary(censor_mod)
cbind(coef(censor_mod), confint(censor_mod))

# summary(w)
# #plot inverse probability weights
# graphics.off()
# ipwplot(weights = w$ipw.weights, timevar = dat4_new$visit,
#         binwidth = 1, ylim = c(-1.5, 1.5), main = "Stabilized inverse probability weights")
# 
# summary(glm(y ~ a, data = dat4_new, weights = w$ipw.weights))

#calculating weights manually
# can also do this with ipwpoint()
# IPT weights:

tmodel1 <- glm(a_1 ~ 1, family=binomial(link=logit), data = dat4)
tmodel2 <- glm(a_2 ~ a_1, family=binomial(link=logit), data = dat4)
# tmodel2_package <- glm(a_2 ~ 1, family=binomial(link=logit), data = dat4)
cmodel <- glm(c ~ w1 + w2 + L1_1 + L2_1 + a_1, family=binomial(link=logit), data = dat4)


smodel1 <- glm(a_1 ~ w1 + w2 + L1_1 + L2_1,
               family=binomial(link=logit), data = dat4)
smodel2 <- glm(a_2 ~  w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1,
               family=binomial(link=logit), data = dat4) #kitchen sink approach
# smodel2_incorrect <- glm(a_2 ~  w1 + w2 +  L1_2 + L2_2,
#                          family=binomial(link=logit), data = dat4)
csmodel <- glm(c ~ a_1,
               family=binomial(link=logit), data = dat4)

# to be updated

# tlp1 <- as.matrix(cbind(1.0)) %*% as.matrix(coef(tmodel1))
# tlp2 <- as.matrix(cbind(1.0, dat4$a_1)) %*% as.matrix(coef(tmodel2))
# 
# slp1 <- as.matrix(cbind(1.0, dat4$w1, dat4$w2, dat4$L1_1, dat4$L2_1)) %*% 
#   as.matrix(coef(smodel1))
# slp2 <- as.matrix(cbind(1.0, dat4$w1, dat4$w2, dat4$L1_1, dat4$L2_1, 
#                         dat4$L1_2, dat4$L2_2, dat4$a_1)) %*% 
#   as.matrix(coef(smodel2))
# 
# num <- tlp2 %*% tlp1 # has negative weights
# deno <- slp1 * slp2

num_test <- predict(tmodel1, type = "response", newdata = dat4)*
  predict(tmodel2, type = "response", newdata = dat4)

deno_test <- predict(smodel1, type = "response", newdata = dat4)*
  predict(smodel2, type = "response", newdata = dat4)

correct_weights <- num_test/deno_test

# fit separate model to combine with weight
cmodel <- glm(c ~ w1 + w2 + L1_1 + L2_1 + a_1, family=binomial(link=logit), data = dat4)
csmodel <- glm(c ~ a_1, family=binomial(link=logit), data = dat4)

cw <- predict(csmodel, type = "response", newdata = dat4)/predict(cmodel, type = "response", newdata = dat4)
w <- correct_weights * cw

censor_design <- svydesign(id=~1, weights = ~ w, data = dat4)
censor_mod <- svyglm(y ~ as.factor(cum_a), design = censor_design)
summary(censor_mod)
cbind(coef(censor_mod), confint(censor_mod))

# num_package <-predict(tmodel1, type = "response", data = dat4)*
#   predict(tmodel2_package, type = "response", data = dat4)
# 
# deno_package<- predict(smodel1, type = "response", data = dat4)*
#   predict(smodel2_incorrect, type = "response", data = dat4)

# weights_package <- num_package/deno_package

# correct_weights <- num_test/deno_test
# incorrect_weights <- num/deno_incorrect

# ds <- svydesign(id=~1, weights = ~ correct_weights, data = dat4)
# summary(svyglm(y ~ as.factor(cum_a), design=ds))

# # calculate using the package's process
# 
# s <- split(w$ipw.weights, 1:2)
# v1 <- as.vector(s[[1]])
# v2 <- as.vector(s[[2]])
# m <- cbind(v1, v2)
# 
# ds_package <- svydesign(id=~1, weights = ~ m[,2], data = dat4)
# summary(svyglm(y ~ as.factor(cum_a), design=ds_package))
# 
# comparison <- tibble(cbind(m,  weights_package, correct_weights))
# View(comparison)

# Weightit package

# examine the initial imbalance at each time point and overall:
bal.tab(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
             a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
        data = dat4, stats = c("m", "ks"), thresholds = c(m = .05),
        which.time = .all)

# specify our weight models for treatment
Wmsm.out <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                             a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                        data = dat4, method = "ps",
                        stabilize = TRUE)
# Error: No missing values are allowed in the treatment variable. Missing values found in a_2.

# take a look at the quality of the weights
summary(Wmsm.out)

# censoring

Wcmsm.out <- weightitMSM(list(c ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                        data = dat4, method = "ps",
                        stabilize = TRUE)

# take a look at the quality of the weights
summary(Wcmsm.out)

wc <- Wcmsm.out$weights * Wmsm.out$weights

# examine how well they perform with respect to covariate balance
bal.tab(Wmsm.out, stats = c("m", "ks"), thresholds = c(m = .05),
        which.time = .none)

# estimate treatment effect
d.w.msm <- svydesign(~1, weights = wc,
                     data = dat4)

cum.fit <- svyglm(y ~ as.factor(cum_a), design = d.w.msm)
summary(cum.fit)

# Confoundr package

# DEMONSTRATION OF DIAGNOSTIC 1

# require censor indicator to be a binary variable and has root name
# require time indicator in the censor name, and must be the same as covariates and exposure
# i.e. if there are two visits, there should be two censoring indicators

dat4 <- dat4 %>% 
  mutate(c_1 = 0,
         c_2 = case_when(c == "censored" ~ 1, TRUE ~ 0)) %>% 
  select(-c)

# create the appropriate history strata
dat4.history <- makehistory.one(
  input=dat4,
  id="id",
  exposure="a",
  times=c(1,2),
  name.history="h"
)

# create a tidy dataframe (long data)
dat4.tidy <- lengthen(
  input=dat4.history,
  id="id",
  diagnostic=1,
  censoring="yes",
  censor = "c",
  exposure="a",
  temporal.covariate=c("L1", "L2"),
  static.covariate = c("w1", "w2"),
  times.exposure=c(1,2), 
  times.covariate=c(1,2),
  history="h"
)

# create a covariate balance table
baltbl4 <- balance(
  input=dat4.tidy,
  diagnostic=1,
  approach="none",
  censoring="yes",
  scope="recent",
  recency=0,
  exposure="a",
  history="h",
  times.exposure=c(1,2),
  times.covariate=c(1,2),
  sd.ref="yes"
)

# create a trellised covariate balance plot
balplot4 <- makeplot(
  input=baltbl4,
  diagnostic=1,
  approach="none",
  censoring="yes",
  scope="recent",
  label.exposure="Treatment",
  label.covariate="Covariate"
)
balplot4

# DEMONSTRATION OF DIAGNOSTIC 3

dat4.history$w_1 <- weights_v1$ipw.weights
dat4.history$w_2 <- weights_v2$ipw.weights

dat4.history3 <- confoundr::lengthen(
  input=dat4.history,
  id="id",
  diagnostic=3,
  censoring="yes",
  censor = "c",
  exposure="a",
  temporal.covariate=c("L1", "L2"),
  static.covariate = c("w1", "w2"),
  times.exposure=c(1,2), 
  times.covariate=c(1,2),
  history="h",
  weight.exposure="w"
)

baltbl4_diag3 <- confoundr::balance(
  input=dat4.history3,
  diagnostic=3,
  approach="weight",
  censoring="yes",
  scope="all",
  exposure="a",
  history="h",
  times.exposure=c(1,2),
  times.covariate=c(1,2),
  sd.ref="yes",
  weight.exposure="w"
)

balplot4_diag3 <- confoundr::makeplot(
  input=baltbl4_diag3,
  diagnostic=3,
  approach="weight",
  censoring="yes",
  scope="all",
  label.exposure="Study Dropout",
  label.covariate="Covariate"
)

balplot4_diag3

# # End

