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

# check this
to_prob <- function(x){
  x <- exp(x)
  x <- x/(x+1)
  return(x)
}
 
y <- rbinom(ntot, 1, prob = to_prob(- 1 + log(0.9)*a_1 + log(0.85)*a_2 + log(1.01)*a_1*a_2 
                            - 0.2*w1 + 0.1*w2 + 0.4*L1_2 + 0.2*L2_2))
# look at simulating binary from odds ratios, different ways to do it
table(y);

# might be worth to try scenarios with larger confounder effect, e.g., - 0.2*a_1 - 0.5*a_2 + 0.01*a_1*a_2;

#2.99 conditional treatment effect of 11 vs 00;

id <- rep(1:1000)
dat2 <- data.frame(id, w1, w2, L1_1, L2_1, a_1, L1_2, L2_2, a_2, y)
# dat2 <- dat2 %>% mutate(cum_a = a_1 + a_2) %>% select(id, cum_a, everything())

# different outcome!
# dat2 <- data.frame(w1, w2, L1_1, L2_1, a_1, L1_2, L2_2, a_2, y)
# length(dat2$y[dat2$y==1])/2000 # proportion = 0.351

# type = "link" for log odds for each person

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
  # calc odds of y=1
  truemean = to_prob(- 1 + log(0.9)*a_1s + log(0.85)*a_2s + log(1.01)*a_1s*a_2s 
                     - 0.2*w1sim + 0.1*w2sim + 0.4*L1_2sim + 0.2*L2_2sim)
  return(mean(truemean))
}

# convert prob into odds

mytrue(a_1s = 1, a_2s = 1) #potential outcome Y(a_1=1, a_2=1);
mytrue(a_1s = 0, a_2s = 0)
mytrue(a_1s = 1, a_2s = 0)
mytrue(a_1s = 0, a_2s = 1)
mytrue(a_1s = 1, a_2s = 1)/mytrue(a_1s = 0, a_2s = 0) # odds ratio

library(arm)
library(ltmle)
library(ipw) #MSM;
library(gfoRmula) #g-computation;

# tmle without superlearner + kitchen sink gform + null Qform (kitchen sink);
tmle_model <- ltmle(dat2,
                    Anodes = c ("a_1","a_2") ,
                    Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                    Ynodes = c("y"), 
                    survivalOutcome =FALSE,
                    # Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a_1",
                    #            y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a_1 + a_2 + a_1*a_2"),
                    gform = c("a_1 ~ w1 + w2 + L1_1 + L2_1",
                              "a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1"),
                    # gcomp = TRUE,
                    # iptw.only = FALSE,
                    # variance.method = "tmle",
                    # SL.library = "default", #with superlearner;
                    abar = list(c(1,1), c(0,0)),
                    estimate.time = FALSE)

summary(tmle_model, estimator="tmle")

# tmle with superlearner + kitchen sink gform + null Qform;
tmle_model_sl <- ltmle(dat2,
                    Anodes = c ("a_1","a_2") ,
                    Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                    Ynodes = c("y"), 
                    survivalOutcome =FALSE,
                    # Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a_1",
                    #            y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a_1 + a_2 + a_1*a_2"),
                    gform = c("a_1 ~ w1 + w2 + L1_1 + L2_1",
                              "a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1"),
                    # gcomp = TRUE,
                    # iptw.only = FALSE,
                    # variance.method = "tmle",
                    SL.library = "SL.glm", #with superlearner for continuous, SL.glm for binary, SL.random_forest, SL.gbm 
                    abar = list(c(1,1), c(0,0)),
                    estimate.time = FALSE)

summary(tmle_model, estimator="tmle")

#tmle package for gcomputation + kitchen sink gform + null Qform;
tmle_model <- ltmle(dat2,
                    Anodes = c ("a_1","a_2") ,
                    Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                    Ynodes = c("y"), 
                    survivalOutcome =FALSE,
                    # Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a_1",
                    #            y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a_1 + a_2 + a_1*a_2"),
                    gform = c("a_1 ~ w1 + w2 + L1_1 + L2_1",
                              "a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1"),
                    gcomp = TRUE,
                    # iptw.only = FALSE,
                    # variance.method = "tmle",
                    abar = list(c(1,1), c(0,0)),
                    estimate.time = FALSE)

summary(tmle_model)

#tmle package for gcomputation + kitchen sink gform + correct Qform;
tmle_model <- ltmle(dat2,
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


