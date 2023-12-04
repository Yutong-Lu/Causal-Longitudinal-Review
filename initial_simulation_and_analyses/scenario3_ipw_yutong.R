# Scenario 3;
# 1000 patients and 3 visits
# y, an end-of-study binary outcome;
# z, a binary treatment assignment;
# w1 and w2 are two baseline covariates (one continuous, one binary) i.e., age, sex;
# L1 and L2 are two time-dependent covariates (one continuous, one binary);
# censoring = no censoring
# Interaction = interaction treatment effect on the outcome;
# non-linearity = non-linear effect for the outcome model

set.seed(123)
ntot=1000; #500 or 200 as small sample example;
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
y <- rnorm(ntot, mean = (1 - 1*a_1 - 2*a_2 + 0.01*a_1*a_2 - 0.2*w1 + 0.1*w2 + 0.4*L1_2 + 0.2*L2_2 + 
                           exp(L1_2) + 0.02*L2_2^3) , sd = 1)
mean(y);

id <- rep(1:1000)

dat3 <- data.frame(w1, w2, L1_1, L2_1, a_1, L1_2, L2_2, a_2, y, id)

dat3 <- dat3 %>% mutate(cum_a = a_1 + a_2)

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
  truemean = (1 - 1*a_1s - 2*a_2s + 0.01*a_1s*a_2s - 0.2*w1sim + 0.1*w2sim + 0.4*L1_2sim + 0.2*L2_2sim + 
                exp(L1_2sim) + 0.02*L2_2sim^3)
  return(mean(truemean))
}

mytrue(a_1s = 1, a_2s = 1) #potential outcome Y(a_1=1, a_2=1);
mytrue(a_1s = 0, a_2s = 0)
mytrue(a_1s = 1, a_2s = 0)
mytrue(a_1s = 0, a_2s = 1)
mytrue(a_1s = 1, a_2s = 1) - mytrue(a_1s = 0, a_2s = 0)

library(tidyverse)
library(arm)
library(ltmle)
library(ipw) #MSM;
library(gfoRmula)

# ipw package

# convert one row per subject to two rows
dat3_new <- dat3 %>%
  pivot_longer(cols = -c(w1,w2,y), 
               names_to = c("variable","visit"), 
               names_sep = "_", 
               values_to = "value") %>% 
  pivot_wider(names_from = variable, values_from = value) %>% 
  mutate(patient = rep(1:1000, each = 2))

dat3_new$visit <- as.numeric(dat3_new$visit)
dat3_new$visit[dat3_new$visit == 1] <- 0
dat3_new$visit[dat3_new$visit == 2] <- 1
dat3_new$y[dat3_new$visit == 0] <- NA
dat3_new$w1 <- as.numeric(dat3_new$w1)
dat3_new$patient <- as.numeric(dat3_new$patient)

glimpse(dat3_new)

w <- ipwtm(
  exposure = a,
  family = "binomial",
  link = "logit",
  numerator = ~ w1 + w2,
  denominator = ~ w1 + w2 + L1 + L2,
  id = patient,
  timevar = visit,
  type = "all",
  data = as.data.frame(dat3_new))

summary(w)
#plot inverse probability weights
graphics.off()
ipwplot(weights = w$ipw.weights, timevar = dat3_new$visit,
        binwidth = 1, ylim = c(-1.5, 1.5), main = "Stabilized inverse probability weights")

summary(glm(y ~ a, data = dat3_new, weights = w$ipw.weights))

# not consistent with simulation!

# manual weight calculation

tmodel1 <- glm(a_1 ~ 1, family=binomial(link=logit), data = dat3)
tmodel2 <- glm(a_2 ~ a_1, family=binomial(link=logit), data = dat3)

smodel1 <- glm(a_1 ~ w1 + w2 + L1_1 + L2_1, 
               family=binomial(link=logit), data = dat3)
smodel2 <- glm(a_2 ~  w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1, 
               family=binomial(link=logit), data = dat3)
smodel2_incorrect <- glm(a_2 ~  w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2, 
                         family=binomial(link=logit), data = dat3)

num <- predict(tmodel1, type = "response", data = dat3)* 
  predict(tmodel2, type = "response", data = dat3)
deno <- predict(smodel1, type = "response", data = dat3)* 
  predict(smodel2, type = "response", data = dat3)
deno_incorrect <- predict(smodel1, type = "response", data = dat3)* 
  predict(smodel2_incorrect, type = "response", data = dat3)

correct_weights <- num/deno
incorrect_weights <- num/deno_incorrect

library(survey)
ds <- svydesign(id=~1, weights = ~ correct_weights, data = dat3)
summary(svyglm(y ~ as.factor(cum_a), design=ds))

ds_incorrect <- svydesign(id=~1, weights = ~ incorrect_weights, data = dat3)
summary(svyglm(y ~ as.factor(cum_a), design=ds_incorrect))



