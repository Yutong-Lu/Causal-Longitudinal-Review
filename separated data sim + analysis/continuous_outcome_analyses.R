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
library(cobalt)
library(confoundr)
options(scipen=999)

install.packages("cobalt")

dat1 <- read_csv("datasets/continuous_outcome_data.csv")


# gfoRmula package;

# creating long format data;
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

summary(gform_cont_eof)

# ltmle packages;

# ltmle without superlearner + kitchen sink gform + Qform (kitchen sink);
tmle_model <- ltmle(dat1,
                    Anodes = c ("a_1","a_2") ,
                    Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                    Ynodes = c("y"), 
                    survivalOutcome =FALSE,
                    Qform = c(L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a_1",
                              y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a_1 + a_2 + a_1*a_2"),
                    gform = c("a_1 ~ w1 + w2 + L1_1 + L2_1",
                              "a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1"),
                    abar = list(c(1,1), c(0,0)),
                    estimate.time = FALSE)

summary(tmle_model, estimator="tmle")

# ltmle with superlearner + kitchen sink gform + Qform;
set.seed(123)
tmle_model <- ltmle(dat1,
                    Anodes = c ("a_1","a_2") ,
                    Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                    Ynodes = c("y"), 
                    survivalOutcome =FALSE,
                    Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a_1",
                               # L2_2 = "Q.kplus1 ~ w1 + w2 + L2_1 + a_1",
                               y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a_1 + a_2 + a_1*a_2"),
                    gform = c("a_1 ~ w1 + w2 + L1_1 + L2_1",
                              "a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1"),
                    SL.library = "default", #with superlearner;
                    abar = list(c(1,1), c(0,0)),
                    estimate.time = FALSE)

summary(tmle_model, estimator="tmle")

# ltmle package for gcomputation + kitchen sink gform + correct Qform;
tmle_model <- ltmle(dat1,
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

summary(tmle_model)

# ipw package;

# creating long format data;
dat1 <- dat1 %>% 
  mutate(id = rep(1:1000),
         cum_a = a_1 + a_2)

dat1_new <- dat1 %>%
  pivot_longer(cols = -c(w1,w2,y,id,cum_a), 
               names_to = c("variable","visit"), 
               names_sep = "_", 
               values_to = "value") %>% 
  pivot_wider(names_from = variable, values_from = value) %>% 
  mutate(time = case_when(visit == 1 ~ 0,
                          visit == 2 ~ 1))

dat1_new$y[dat1_new$visit == 1] <- NA

dat1_new$visit <- as.numeric(dat1_new$visit)
dat1_new$w1 <- as.numeric(dat1_new$w1)

dat1_new$cum_a[dat1_new$time == 0] <- dat1$a_1
dat1_new$cum_a[dat1_new$time == 1] <- dat1$a_1+dat1$a_2

w <- ipwtm(
  exposure = a,
  family = "binomial",
  link = "logit",
  numerator = ~ w1+ w2,
  denominator = ~ w1 + w2 + L1 + L2,
  id = id,
  timevar = visit,
  type = "all",
  data = as.data.frame(dat1_new))

# calculate using the package's process

s <- split(w$ipw.weights, 1:2)
v1 <- as.vector(s[[1]])
v2 <- as.vector(s[[2]])
m <- cbind(v1, v2)

ds_package <- svydesign(id=~1, weights = ~ m[,2], data = dat1)
summary(svyglm(y ~ as.factor(cum_a), design=ds_package))
confint(svyglm(y ~ as.factor(cum_a), design=ds_package))
# simulation -> no huge confounding, no huge bias
# true empirical studies many simulation

summary(w)
#plot inverse probability weights
graphics.off()
ipwplot(weights = w$ipw.weights, timevar = dat1_new$visit,
        binwidth = 1, ylim = c(-1.5, 1.5), 
        main = "Stabilized inverse probability weights")

summary(glm(y ~ a, data = dat1_new, weights = w$ipw.weights))

# model using ipwpoint

weights_v1 <- ipwpoint(
  exposure = a_1,
  family = "gaussian",
  numerator = ~ 1,
  denominator = ~ w1 + w2 + L1_1 + L2_1,
  data = as.data.frame(dat1))

weights_v2 <- ipwpoint(
  exposure = a_2,
  family = "gaussian",
  numerator = ~ a_1,
  denominator = ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1,
  data = as.data.frame(dat1))

w <- weights_v1$ipw.weights * weights_v2$ipw.weights

cont_design <- svydesign(id=~1, weights = ~ w, data = dat1)
cont_mod <- svyglm(y ~ as.factor(cum_a), design = cont_design)
summary(cont_mod)
cbind(coef(cont_mod), confint(cont_mod))

# #calculating weights manually
# # IPT weights:
# 
# tmodel1 <- glm(a_1 ~ 1, family=binomial(link=logit), data = dat1)
# tmodel2 <- glm(a_2 ~ a_1, family=binomial(link=logit), data = dat1)
# tmodel2_package <- glm(a_2 ~ 1, family=binomial(link=logit), data = dat1)
# 
# 
# smodel1 <- glm(a_1 ~ w1 + w2 + L1_1 + L2_1, 
#                family=binomial(link=logit), data = dat1)
# smodel2 <- glm(a_2 ~  w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1, 
#                family=binomial(link=logit), data = dat1) #kitchen sink approach
# smodel2_incorrect <- glm(a_2 ~  w1 + w2 +  L1_2 + L2_2, 
#                          family=binomial(link=logit), data = dat1)
# 
# num <- predict(tmodel1, type = "response", data = dat1)* 
#   predict(tmodel2, type = "response", data = dat1)
# deno <- predict(smodel1, type = "response", data = dat1)* 
#   predict(smodel2, type = "response", data = dat1)
# deno_incorrect <- predict(smodel1, type = "response", data = dat1)* 
#   predict(smodel2_incorrect, type = "response", data = dat1)
# 
# num_package <-predict(tmodel1, type = "response", data = dat1)* 
#   predict(tmodel2_package, type = "response", data = dat1)
# 
# deno_package<- predict(smodel1, type = "response", data = dat1)* 
#   predict(smodel2_incorrect, type = "response", data = dat1)
# 
# weights_package <- num_package/deno_package
# 
# correct_weights <- num/deno
# incorrect_weights <- num/deno_incorrect
# 
# ds <- svydesign(id=~1, weights = ~ correct_weights, data = dat1)
# summary(svyglm(y ~ as.factor(cum_a), design=ds))
# confint(svyglm(y ~ as.factor(cum_a), design=ds))
# 
# 
# comparison <- tibble(cbind(m,  weights_package, correct_weights))
# View(comparison)

# Weightit package

# examine the initial imbalance at each time point and overall:
bal.tab(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
             a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
        data = dat1, stats = c("m", "ks"), thresholds = c(m = .05),
        which.time = .all)

# specify our weight models
Wmsm.out <- weightitMSM(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                        a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                        data = dat1, method = "ps",
                        stabilize = TRUE)

# take a look at the quality of the weights 
summary(Wmsm.out)

# examine how well they perform with respect to covariate balance
bal.tab(Wmsm.out, stats = c("m", "ks"), thresholds = c(m = .05),
        which.time = .none)

# estimate treatment effect
d.w.msm <- svydesign(~1, weights = Wmsm.out$weights,
                     data = dat1)

cum.fit <- svyglm(y ~ as.factor(cum_a), design = d.w.msm)
summary(cum.fit)

# Confoundr package

# DEMONSTRATION OF DIAGNOSTIC 1

# create the appropriate history strata
dat1.history <- makehistory.one(
  input=dat1,
  id="id",
  exposure="a",
  times=c(1,2),
  name.history="h"
)

# create a tidy dataframe (long data)
dat1.tidy <- confoundr::lengthen(
  input=dat1.history,
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
baltbl1 <- confoundr::balance(
  input=dat1.tidy,
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
balplot1 <- confoundr::makeplot(
  input=baltbl1,
  diagnostic=1,
  approach="none",
  censoring="no",
  scope="recent",
  label.exposure="Treatment",
  label.covariate="Covariate"
)
balplot1

# DEMONSTRATION OF DIAGNOSTIC 3

dat1.history$w_1 <- weights_v1$ipw.weights
dat1.history$w_2 <- weights_v2$ipw.weights

dat1.history3 <- confoundr::lengthen(
  input=dat1.history,
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

baltbl1_diag3 <- confoundr::balance(
  input=dat1.history3,
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

balplot1_diag3 <- makeplot(
  input=baltbl1_diag3,
  diagnostic=3,
  approach="weight",
  censoring="no",
  scope="all",
  label.exposure="Study Dropout",
  label.covariate="Covariate"
)

balplot1_diag3

# # End
