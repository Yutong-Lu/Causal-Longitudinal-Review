# Scenario 1;
# 1000 patients and 3 visits
# y, an end-of-study continuous outcome;
# z, a binary treatment assignment;
# w1 and w2 are two baseline covariates (one continuous, one binary) i.e., age, sex;
# L1 and L2 are two time-dependent covariates (one continuous, one binary);
# censoring = no censoring
# Interaction = interaction treatment effect on the outcome;
# linearity = additive linear effect for the outcome model;

set.seed(123)
ntot=1000;
samplesize=10000

# visit 1;
# covariates;
w1 <- rbinom(ntot, 1, prob = 0.5) #50-50 male and female (reference);
w2 <- rnorm(ntot, mean = 12, sd = 4) #age;
L1_1 <- rbinom(ntot, 1, prob = plogis(-1 + 0.01*w1 + 0.01*w2)) #a age & sex related baseline binary clinical variable;
L2_1 <- rnorm(ntot, mean = (0.01*w1 + 0.01*w2), sd = 1) #a age & sex related baseline continuous clinical variable;

# Exposure;
a_1 <- rbinom(ntot, 1, plogis(-1 - 0.01*w1 + 0.01*w2 - 0.1*L1_1 + 0.1*L2_1))
table(a_1)

# visit 2;
L1_2 <- rbinom(ntot, 1, plogis(-1 + 0.01*w1 + 0.01*w2 + 0.2*L1_1 - 0.2*a_1))
L2_2 <- rnorm(ntot, mean = (L2_1 - 0.01*w1 + 0.01*w2 - 0.2*a_1), sd = 1)
a_2 <- rbinom(ntot, 1, plogis(-1 - 0.01*w1 + 0.01*w2 - 0.1*L1_2 + 0.1*L2_2 + a_1))
table(a_1, a_2)

# visit 3;
# disease activity score, continuous;
# a_1 and a_2 are dominant predictors for y, approximately 1 unit increase in a_1 result in 1 unit change on y
# but confounders have smaller effect in this setting, 1 unit increase in w2, y is only changed on average by 0.1;
y <- rnorm(ntot, mean = (1 - 1*a_1 - 2*a_2 + 0.01*a_1*a_2 - 0.2*w1 + 0.1*w2 + 0.4*L1_2 + 0.2*L2_2) , sd = 1)
mean(y);

# might be worth to try scenarios with larger confounder effect, e.g., - 0.2*a_1 - 0.5*a_2 + 0.01*a_1*a_2;

#2.99 conditional treatment effect of 11 vs 00;

id <- rep(1:1000)

dat1 <- data.frame(w1, w2, L1_1, L2_1, a_1, L1_2, L2_2, a_2, y, id)

dat1 <- dat1 %>% mutate(cum_a = a_1 + a_2)

long <- data.frame(w1=c(w1, w1), w2=c(w2,w2), L1 = c(L1_1, L1_2), L2 = c(L2_1, L2_2),
                   y= c(rep(NA, 1000), y), a = c(a_1, a_2), cum_a = c(a_1, a_1+a_2), 
                   visit=c(rep(1,1000), rep(2,1000)),
                   id = c(rep(1:1000,2)))

mytrue <- function(a_1s = 1, a_2s = 1){
  #visit 1;
  w1sim <- rbinom(samplesize, 1, prob = 0.5) #50-50 male and female;
  w2sim <- rnorm(samplesize, mean = 12, sd = 4) #age;
  L1_1sim <- rbinom(samplesize, 1, prob = plogis(-1 + 0.01*w1sim + 0.01*w2sim)) #a age related baseline binary clinical variable;
  L2_1sim <- rnorm(samplesize, mean = (0.01*w1sim + 0.01*w2sim), sd = 1) #a age & sex related baseline continuous clinical variable;
  
  #Visit 2, simulate all potential variables;
  L1_2sim <- rbinom(samplesize, 1, plogis(-1 + 0.01*w1sim + 0.01*w2sim + 0.2*L1_1sim - 0.2*a_1s))
  L2_2sim <- rnorm(samplesize, mean = (L2_1sim - 0.01*w1sim + 0.01*w2sim - 0.2*a_1s), sd = 1)
  
  #Visit 3, simulate potential outcomes;
  truemean = (1 - 1*a_1s - 2*a_2s + 0.01*a_1s*a_2s - 0.1*w1sim + 0.01*w2sim + 0.1*L1_2sim + 0.1*L2_2sim)
  return(mean(truemean))
}

mytrue(a_1s = 1, a_2s = 1) #potential outcome Y(a_1=1, a_2=1);
mytrue(a_1s = 0, a_2s = 0)
mytrue(a_1s = 1, a_2s = 0)
mytrue(a_1s = 0, a_2s = 1)
mytrue(a_1s = 1, a_2s = 1) - mytrue(a_1s = 0, a_2s = 0); #-3.010579, this can change due to random sampling;

library(tidyverse)
library(arm)
library(ltmle)
library(ipw) #MSM;
library(gfoRmula) #g-computation;

# ipw package

# convert one row per subject to two rows
dat_new <- dat1 %>%
  pivot_longer(cols = -c(w1,w2,y,id), 
               names_to = c("variable","visit"), 
               names_sep = "_", 
               values_to = "value") %>% 
  pivot_wider(names_from = variable, values_from = value) %>% 
  mutate(patient = rep(1:1000, each = 2))

dat_new$visit <- as.numeric(dat_new$visit)
dat_new$time[dat_new$visit == 1] <- 0
dat_new$time[dat_new$visit == 2] <- 1
dat_new$y[dat_new$visit == 1] <- NA
dat_new$w1 <- as.numeric(dat_new$w1)
dat_new$patient <- as.numeric(dat_new$patient)

# create a cumulative a to that time point
# v1 cum_a <- a_1
# v2 cum_a <- sum(a_1+a_2)

dat_new$cum_a[dat_new$time == 0] <- dat1$a_1
dat_new$cum_a[dat_new$time == 1] <- dat1$a_1+dat1$a_2

glimpse(dat_new)

# w <- ipwtm(
#   exposure = c(a_1, a_2),
#   family = "binomial",
#   link = "logit",
#   numerator = c("~ w1 + w2", "~ w1 + w2 + a_1"),
#   denominator = c("~ w1 + w2 + L1_1 + L1_2 + L2_1 + L2_2", "~ w1 + w2 + L1_1 + L1_2 
#                   + L2_1 + L2_2 + a_1"),
#   id = patient,
#   timevar = visit,
#   type = "all",
#   data = as.data.frame(dat1))

w <- ipwtm(
  exposure = a,
  family = "binomial",
  link = "logit",
  numerator = ~ w1+ w2,
  denominator = ~ w1 + w2 + L1 + L2,
  id = patient,
  timevar = visit,
  type = "all",
  data = as.data.frame(dat_new))

summary(w)
#plot inverse probability weights
graphics.off()
ipwplot(weights = w$ipw.weights, timevar = dat_new$visit,
        binwidth = 1, ylim = c(-1.5, 1.5), 
        main = "Stabilized inverse probability weights")

summary(glm(y ~ a, data = dat_new, weights = w$ipw.weights))

#calculating weights manually
# IPT weights:

tmodel1 <- glm(a_1 ~ 1, family=binomial(link=logit), data = dat1)
tmodel2 <- glm(a_2 ~ a_1, family=binomial(link=logit), data = dat1)

# tlp1 <- as.matrix(cbind(1.0, obs[,1:2])) %*% as.matrix(coef(tmodel1))
# tlp2 <- as.matrix(cbind(1.0, obs[,3:4], z[,1])) %*% as.matrix(coef(tmodel2))

smodel1 <- glm(a_1 ~ w1 + w2 + L1_1 + L2_1, 
               family=binomial(link=logit), data = dat1)
smodel2 <- glm(a_2 ~  w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1, 
               family=binomial(link=logit), data = dat1)
smodel2_incorrect <- glm(a_2 ~  w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2, 
               family=binomial(link=logit), data = dat1)

# slp1 <- as.matrix(cbind(rep(1.0, ntot))) %*% as.matrix(coef(smodel1))
# slp2 <- as.matrix(cbind(1.0, z[,1])) %*% as.matrix(coef(smodel2))

# pt <- (exp(z[,1] * tlp1)/(1+exp(tlp1))) * (exp(z[,2] * tlp2)/(1+exp(tlp2)))
# 
# sc <- (exp(z[,1] * slp1)/(1+exp(slp1))) * (exp(z[,2] * slp2)/(1+exp(slp2)))

num <- predict(tmodel1, type = "response", data = dat1)* 
  predict(tmodel2, type = "response", data = dat1)
deno <- predict(smodel1, type = "response", data = dat1)* 
  predict(smodel2, type = "response", data = dat1)
deno_incorrect <- predict(smodel1, type = "response", data = dat1)* 
  predict(smodel2_incorrect, type = "response", data = dat1)

# iptw1 <- 1.0/(exp(z[,1] * tlp1)/(1+exp(tlp1)))
# iptws1 <- (exp(z[,1] * slp1)/(1+exp(slp1))) * iptw1
# 
# iptw2 <- 1.0/pt
# iptws2 <- sc * iptw2

correct_weights <- num/deno
incorrect_weights <- num/deno_incorrect

library(survey)
ds <- svydesign(id=~1, weights = ~ correct_weights, data = dat1)
summary(svyglm(y ~ as.factor(cum_a), design=ds))

ds_incorrect <- svydesign(id=~1, weights = ~ incorrect_weights, data = dat1)
summary(svyglm(y ~ as.factor(cum_a), design=ds_incorrect))


# create a matrix of n rows corresponding patient ids (1000 rows)
# col: 2 visits (weight for individual i at visit 1 and 2)
# ipw.weight format?
# odd index col 1 and even index col 2
# then use gee 

# library(geepack)
# 
# geese(y ~ as.factor(cum_a), id = id, family = "gaussian", data = long, 
#       corstr = "independence", weights = v2)
# weights doesnt work here

# w$ipw.weights
# s <- split(w$ipw.weights, 1:2)
# v1 <- as.vector(s[[1]])
# v2 <- as.vector(s[[2]])
# m <- cbind(v1, v2)
# m

#doesnt allow separate visit specification, require a long dataset

comparison <- tibble(cbind(m, correct_weights, incorrect_weights))
