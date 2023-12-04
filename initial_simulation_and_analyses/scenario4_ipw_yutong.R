# Scenario 4;
# 400 patients and 3 visits
# y, an end-of-study continuous outcome;
# z, a binary treatment assignment;
# w1 and w2 are two baseline covariates (one continuous, one binary) i.e., age, sex;
# L1 and L2 are two time-dependent covariates (one continuous, one binary);
# censoring = has censoring
# Interaction = interaction treatment effect on the outcome;
# linearity = additive linear effect for the outcome model;

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
y <- rnorm(ntot, mean = (1 - 1*a_1 - 2*a_2 + 0.01*a_1*a_2 - 0.2*w1 + 0.1*w2 + 0.4*L1_2 + 0.2*L2_2) , sd = 1)
mean(y);

# censoring
# simulate c from a binary+logist format;
c <- rbinom(ntot, 1, plogis(-1 + 0.02*w1 + 0.01*w2 - 0.4*L1_1 - 0.2*L2_1 - 0.5*a_1))

dat4 <- data.frame(w1, w2, L1_1, L2_1, a_1, c, L1_2, L2_2, a_2, y)
sum(dat4$c==1)/2000 #0.122

dat4$L1_2[dat4$c==1] <- NA
dat4$L2_2[dat4$c==1] <- NA
dat4$a_2[dat4$c==1] <- NA
dat4$y[dat4$c==1] <- NA

# dat4$c[dat4$c==0] <- "uncensored"
# dat4$c[dat4$c==1] <- "censored"
# dat4$c <- as.factor(dat4$c)

complete_dat4 <- dat4[complete.cases(dat4), ]

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

library(arm)
library(ltmle)
library(ipw) #MSM;
library(gfoRmula) 

# convert one row per subject to two rows
dat4_new <- dat4 %>%
  pivot_longer(cols = -c(w1,w2,y,c), 
               names_to = c("variable","visit"), 
               names_sep = "_", 
               values_to = "value") %>% 
  pivot_wider(names_from = variable, values_from = value) %>% 
  mutate(patient = rep(1:1000, each = 2))

dat4_new$visit <- as.numeric(dat4_new$visit)
dat4_new$visit[dat4_new$visit == 1] <- 0
dat4_new$visit[dat4_new$visit == 2] <- 1
dat4_new$y[dat4_new$visit == 0] <- NA
dat4_new$w1 <- as.numeric(dat4_new$w1)
dat4_new$patient <- as.numeric(dat4_new$patient)

# ipw

w <- ipwtm(
  exposure = a,
  family = "binomial",
  link = "logit",
  numerator = ~ w1 + w2,
  denominator = ~ w1 + w2 + L1 + L2,
  id = patient,
  timevar = visit,
  type = "all",
  data = as.data.frame(dat4_new))

summary(w)
#plot inverse probability weights
graphics.off()
ipwplot(weights = w$ipw.weights, timevar = dat4_new$visit,
        binwidth = 1, ylim = c(-1.5, 1.5), main = "Stabilized inverse probability weights")

summary(glm(y ~ a, data = dat4_new, weights = w$ipw.weights))

# manual weight calculation

tmodel1 <- glm(a_1 ~ 1, family=binomial(link=logit), data = dat4)
tmodel2 <- glm(a_2 ~ a_1, family=binomial(link=logit), data = dat4)

smodel1 <- glm(a_1 ~ w1 + w2 + L1_1 + L2_1, 
               family=binomial(link=logit), data = dat4)
smodel2 <- glm(a_2 ~  w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1, 
               family=binomial(link=logit), data = dat4)
smodel2_incorrect <- glm(a_2 ~  w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2, 
                         family=binomial(link=logit), data = dat4)

num <- predict(tmodel1, type = "response", data = dat4)* 
  predict(tmodel2, type = "response", data = dat4)
deno <- predict(smodel1, type = "response", data = dat4)* 
  predict(smodel2, type = "response", data = dat4)
deno_incorrect <- predict(smodel1, type = "response", data = dat4)* 
  predict(smodel2_incorrect, type = "response", data = dat4)

correct_weights <- num/deno
incorrect_weights <- num/deno_incorrect

library(survey)
ds <- svydesign(id=~1, weights = ~ correct_weights, data = dat4)
summary(svyglm(y ~ as.factor(cum_a), design=ds))

ds_incorrect <- svydesign(id=~1, weights = ~ incorrect_weights, data = dat4)
summary(svyglm(y ~ as.factor(cum_a), design=ds_incorrect))

