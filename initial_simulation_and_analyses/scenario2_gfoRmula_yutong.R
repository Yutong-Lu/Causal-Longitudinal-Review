# Scenario 2;
# 1000 patients and 3 visits
# y, an end-of-study binary outcome;
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
y <- rbinom(ntot, 1, plogis(1 - 1*a_1 - 2*a_2 + 0.01*a_1*a_2 - 0.2*w1 + 0.1*w2 + 0.4*L1_2 + 0.2*L2_2))
mean(y);

# might be worth to try scenarios with larger confounder effect, e.g., - 0.2*a_1 - 0.5*a_2 + 0.01*a_1*a_2;

#2.99 conditional treatment effect of 11 vs 00;

dat2 <- data.frame(w1, w2, L1_1, L2_1, a_1, L1_2, L2_2, a_2, y)
length(dat2$y[dat2$y==1])/2000 # proportion = 0.351

numerator <- sum(dat2$y[dat2$a_1==1&dat2$a_2==1])/sum(dat2$y[dat2$a_1==0&dat2$a_2==0])
denominator <- (length(dat2$y[dat2$a_1==1&dat2$a_2==1])-sum(dat2$y[dat2$a_1==1&dat2$a_2==1]))/
  (length(dat2$y[dat2$a_1==0&dat2$a_2==0])-sum(dat2$y[dat2$a_1==0&dat2$a_2==0]))
or <- numerator/denominator
or #0.05011588

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
  truemean = plogis(1 - 1*a_1s - 2*a_2s + 0.01*a_1s*a_2s - 0.1*w1sim + 0.01*w2sim + 0.1*L1_2sim + 0.1*L2_2sim)
  return(mean(truemean))
}

mytrue(a_1s = 1, a_2s = 1) #potential outcome Y(a_1=1, a_2=1);
mytrue(a_1s = 0, a_2s = 0)
mytrue(a_1s = 1, a_2s = 0)
mytrue(a_1s = 0, a_2s = 1)
mytrue(a_1s = 1, a_2s = 1) - mytrue(a_1s = 0, a_2s = 0) #-0.6211263

library(arm)
library(ltmle)
library(ipw) #MSM;
library(gfoRmula) #g-computation;

# gfoRumla package

dat2_new <- dat2 %>%
  pivot_longer(cols = -c(w1,w2,y), 
               names_to = c("variable","visit"), 
               names_sep = "_", 
               values_to = "value") %>% 
  pivot_wider(names_from = variable, values_from = value) %>% 
  mutate(patient = rep(1:1000,each = 2))

dat2_new$visit <- as.numeric(dat2_new$visit)
dat2_new$visit[dat2_new$visit == 1] <- 0
dat2_new$visit[dat2_new$visit == 2] <- 1
dat2_new$y[dat2_new$visit == 0] <- NA
dat2_new$w1 <- as.numeric(dat2_new$w1)
dat2_new$patient <- as.numeric(dat2_new$patient)

# gfoRmula package

id <- 'patient'
time_name <- 'visit'
covnames <- c("L1", "L2", "a")
outcome_name <- 'y'
covtypes <- c('binary', 'normal', 'binary')
histories <- c(lagged, cumavg)
histvars <- list(c('a', 'L1', 'L2'),c('L1','L2'))

covparams <- list(covmodels = c(L1 ~ w1 + w2 + lag1_L1 + lag1_a,
                                L2 ~ lag1_L2 + w1 + w2 + lag1_a,
                                a ~ w1 + w2 + lag1_L1 + L2 + lag1_a))
ymodel <- y ~ lag1_a + a + lag1_a*a + w1 + w2 + L1 + L2

intvars <- list('a', 'a')
interventions <- list(list(c(static, rep(0, 2))),
                      list(c(threshold, 1, Inf)))
int_descript <- c('Never treat', 'Threshold - lower bound 1')

gform_bin_eof <- gformula_binary_eof(
  obs_data = dat2_new,
  id = id,
  time_name = time_name,
  covnames =covnames,
  outcome_name = outcome_name, 
  covtypes = c("binary", "normal", "binary"),
  covparams = covparams,  ymodel = ymodel,
  intvars = intvars, interventions = interventions,
  int_descript = int_descript, 
  histories = c(lagged), histvars = list(c('a',"L1","L2")),
  basecovs = c("w1","w2"), nsimul = 1000,
  nsamples = 1000, parallel = TRUE, ncores = 2,
  seed = 123)

gform_bin_eof
