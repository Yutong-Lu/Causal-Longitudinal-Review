# Scenario 1;
# 400 patients and 3 visits
# y, an end-of-study continuous outcome;
# z, a binary treatment assignment;
# w1 and w2 are two baseline covariates (one continuous, one binary) i.e., age, sex;
# L1 and L2 are two time-dependent covariates (one continuous, one binary);
# censoring = no censoring
# Interaction = interaction treatment effect on the outcome;
# linearity = additive linear effect for the outcome model;

set.seed(123)
<<<<<<< Updated upstream
ntot=1000; #500 or 200 as small sample example;
samplesize=10000
=======
ntot=1000 # try 500 or 200 as small example
samplesize=10000 
>>>>>>> Stashed changes

# visit 1;
# covariates;
w1 <- rbinom(ntot, 1, prob = 0.5) #50-50 male and female (reference);
w2 <- rnorm(ntot, mean = 12, sd = 4) #age;
# increasing age, increasing prob
L1_1 <- rbinom(ntot, 1, prob = plogis(-1 + 0.01*w1 + 0.01*w2)) #a age & sex related baseline binary clinical variable;
L2_1 <- rnorm(ntot, mean = (0.01*w1 + 0.01*w2), sd = 1) #a age & sex related baseline continuous clinical variable;

# Exposure;
# look at age and sex and current covariates (presented values)
a1 <- rbinom(ntot, 1, plogis(-1 - 0.01*w1 + 0.01*w2 - 0.1*L1_1 + 0.1*L2_1))
table(a1)

# visit 2;
# predicted by the baseline + previous cov and trt
L1_2 <- rbinom(ntot, 1, plogis(-1 + 0.01*w1 + 0.01*w2 + 0.2*L1_1 - 0.2*a1))
L2_2 <- rnorm(ntot, mean = (L2_1 - 0.01*w1 + 0.01*w2 - 0.2*a1), sd = 1)
# switch meds or stay the same
a2 <- rbinom(ntot, 1, plogis(-1 - 0.01*w1 + 0.01*w2 - 0.1*L1_2 + 0.1*L2_2 + a1))
table(a1, a2)

# visit 3;
# disease activity score, continuous;
<<<<<<< Updated upstream
# a1 and a2 are dominant predictors for y, approximately 1 unit increase in a1 result in 1 unit change on y
# but confounders have smaller effect in this setting, 1 unit increase in w2, y is only changed on average by 0.1;
y <- rnorm(ntot, mean = (1 - 1*a1 - 2*a2 + 0.01*a1*a2 - 0.2*w1 + 0.1*w2 + 0.4*L1_2 + 0.2*L2_2) , sd = 1)
mean(y);

# might be worth to try scenarios with larger confounder effect, e.g., - 0.2*a1 - 0.5*a2 + 0.01*a1*a2;

#2.99 conditional treatment effect of 11 vs 00;
=======
# a1 and a2 are dom. pred. for y, approx 1 unit increase in a1 results in 1 unit change on y
# but confounders have smaller effect in this setting, 1 unit increase in w2, y is only changed on avg by 0.01
# create a new simulated data and try to make it a binary (not rare outcome, 
# ideally 200 outcomes, 20%)
# in this example, completely linear relationship, no power to 2
# can consider non-linear effects
y <- rnorm(ntot, mean = (1 - 1*a1 - 2*a2 + 0.01*a1*a2 - 0.1*w1 + 0.01*w2 + 0.1*L1_2 + 0.1*L2_2) , sd = 1)
mean(y)
>>>>>>> Stashed changes

#2.99 condi trt effect of 11 vs 00

# might be worth to try sencarios with larger confounder effects, 
# e.g. -0.2*a1-0.5*a2 + 0.01*a1*a2
# using lower a or higher w

dat1 <- data.frame(w1, w2, L1_1, L2_1, a1, L1_2, L2_2, a2, y)

mytrue <- function(a1s = 1, a2s = 1){
  #visit 1;
  #does not involve treatment assignment, fake population
  w1sim <- rbinom(samplesize, 1, prob = 0.5) #50-50 male and female;
  w2sim <- rnorm(samplesize, mean = 12, sd = 4) #age;
  L1_1sim <- rbinom(samplesize, 1, prob = plogis(-1 + 0.01*w1sim + 0.01*w2sim)) #a age related baseline binary clinical variable;
  L2_1sim <- rnorm(samplesize, mean = (0.01*w1sim + 0.01*w2sim), sd = 1) #a age & sex related baseline continuous clinical variable;

  #Visit 2, simulate all potential variables;
  #baseline + desired trt
  #trick: having everyone in the cohort receive the same trt and doesnt receive
  #causal para: contrast, nothing to do with confounders
  L1_2sim <- rbinom(samplesize, 1, plogis(-1 + 0.01*w1sim + 0.01*w2sim + 0.2*L1_1sim - 0.2*a1s))
  L2_2sim <- rnorm(samplesize, mean = (L2_1sim - 0.01*w1sim + 0.01*w2sim - 0.2*a1s), sd = 1)

  #Visit 3, simulate potential outcomes;
  truemean = (1 - 1*a1s - 2*a2s + 0.01*a1s*a2s - 0.1*w1sim + 0.01*w2sim + 0.1*L1_2sim + 0.1*L2_2sim)
  return(mean(truemean))
}

<<<<<<< Updated upstream
mytrue(a1s = 1, a2s = 1) #potential outcome Y(a1=1, a2=1);
mytrue(a1s = 0, a2s = 0)
mytrue(a1s = 1, a2s = 0)
mytrue(a1s = 0, a2s = 1)
mytrue(a1s = 1, a2s = 1) - mytrue(a1s = 0, a2s = 0); #-3.010579, this can change due to random sampling;

# > mytrue(a1s = 1, a2s = 1) #potential outcome Y(a1=1, a2=1);
# [1] -1.887657
# > mytrue(a1s = 0, a2s = 0)
# [1] 1.125843
# > mytrue(a1s = 1, a2s = 0)
# [1] 0.1039287
# > mytrue(a1s = 0, a2s = 1)
# [1] -0.873694
# > mytrue(a1s = 1, a2s = 1) - mytrue(a1s = 0, a2s = 0);
# [1] -3.010579

library(arm)
library(ltmle)
library(ipw) #MSM;
library(gfoRmula) #g-computation;

# tmle without superlearner + kitchen sink gform + null Qform (kitchen sink);
=======
mytrue(a1s = 1, a2s = 1) #potential outcome y(a1=1,a2=1)
mytrue(a1s = 0, a2s = 0)
mytrue(a1s = 1, a2s = 0)
mytrue(a1s = 0, a2s = 1)
mytrue(a1s = 1, a2s = 1) - mytrue(a1s = 0, a2s = 0) #-3.014203 this can change due to random sampling

library(ltmle)
library(arm) # dependent
library(ipw)
library(gfoRmula)

# tmle without superlearner (on Q) + kitchen sink gform + null Qform;
>>>>>>> Stashed changes
tmle_model <- ltmle(dat1,
                    Anodes = c ("a1","a2") ,
                    Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                    Ynodes = c("y"), 
                    survivalOutcome =FALSE,
                    
                    # default Q is kitchen sink
                    # covaritates mod
                    # Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a1",
                    #            y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a1 + a2 + a1*a2"),
                    
                    # trt assignment func
                    # we didnt consider L11 and L21 in simulation
                    # but we dont know this in the real data -> throw all available variables in
                    # remember there's a time factor: a1 only L1
                    gform = c("a1 ~ w1 + w2 + L1_1 + L2_1",
                              "a2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a1"),
                    
                    # gcomp = TRUE,
                    # iptw.only = FALSE,
                    # variance.method = "tmle",
                    # SL.library = "default", #with superlearner;
                    abar = list(c(1,1), c(0,0)), #potential 11 00 and their contrast
                    estimate.time = FALSE)

# we use glm, not speedglm
# default Q is not the same as the sim truth for both L12 and y

summary(tmle_model, estimator="tmle")
# y00 is a little bit off -> true eff a little bit off but ci coverage is achieved

# tmle with superlearner + kitchen sink gform + null Qform;
# use when there's complex relationships, more robust
# not work well with small sample sizes!
tmle_model <- ltmle(dat1,
                    Anodes = c ("a1","a2") ,
                    Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                    Ynodes = c("y"), 
                    survivalOutcome =FALSE,
                    # Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a1",
                    #            y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a1 + a2 + a1*a2"),
                    gform = c("a1 ~ w1 + w2 + L1_1 + L2_1",
                              "a2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a1"),
                    # gcomp = TRUE,
                    # iptw.only = FALSE,
                    # variance.method = "tmle",
                    SL.library = "default", #with superlearner;
                    abar = list(c(1,1), c(0,0)),
                    estimate.time = FALSE)

summary(tmle_model, estimator="tmle")

#tmle package for gcomputation + kitchen sink gform + null Qform;
tmle_model <- ltmle(dat1,
                    Anodes = c ("a1","a2") ,
                    Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                    Ynodes = c("y"), 
                    survivalOutcome =FALSE,
                    # Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a1",
                    #            y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a1 + a2 + a1*a2"),
                    gform = c("a1 ~ w1 + w2 + L1_1 + L2_1",
                              "a2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a1"),
                    gcomp = TRUE,
                    # iptw.only = FALSE,
                    # variance.method = "tmle",
                    abar = list(c(1,1), c(0,0)),
                    estimate.time = FALSE)

summary(tmle_model)

#tmle package for gcomputation + kitchen sink gform + correct Qform;
tmle_model <- ltmle(dat1,
                    Anodes = c ("a1","a2") ,
                    Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                    Ynodes = c("y"), 
                    survivalOutcome =FALSE,
                    Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a1",
                               y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a1 + a2 + a1*a2"),
                    gform = c("a1 ~ w1 + w2 + L1_1 + L2_1",
                              "a2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a1"),
                    gcomp = TRUE,
                    # iptw.only = FALSE,
                    # variance.method = "tmle",
                    abar = list(c(1,1), c(0,0)),
                    estimate.time = FALSE)

summary(tmle_model)
# better est than null Qform in gcomp

# very bad if ci doesnt cover the true

# figure out how to run ipw and gfoRmula, create a new file
# what simulation senarios to consider

# moderate effect (1/2 of trt)
