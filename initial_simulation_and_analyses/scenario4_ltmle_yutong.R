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
a1 <- rbinom(ntot, 1, plogis(-1 - 0.01*w1 + 0.01*w2 - 0.1*L1_1 + 0.1*L2_1))
table(a1)

# visit 2;
L1_2 <- rbinom(ntot, 1, plogis(-1 + 0.01*w1 + 0.01*w2 + 0.2*L1_1 - 0.2*a1))
L2_2 <- rnorm(ntot, mean = (L2_1 - 0.01*w1 + 0.01*w2 - 0.2*a1), sd = 1)
a2 <- rbinom(ntot, 1, plogis(-1 - 0.01*w1 + 0.01*w2 - 0.1*L1_2 + 0.1*L2_2 + a1))
table(a1, a2)

# visit 3;
# disease activity score, continuous;
# a1 and a2 are dominant predictors for y, approximately 1 unit increase in a1 result in 1 unit change on y
# but confounders have smaller effect in this setting, 1 unit increase in w2, y is only changed on average by 0.1;
y <- rnorm(ntot, mean = (1 - 1*a1 - 2*a2 + 0.01*a1*a2 - 0.2*w1 + 0.1*w2 + 0.4*L1_2 + 0.2*L2_2) , sd = 1)
mean(y);

# censoring
# simulate c from a binary+logist format;
c <- rbinom(ntot, 1, plogis(-1 + 0.02*w1 + 0.01*w2 - 0.4*L1_1 - 0.2*L2_1 - 0.5*a1))

dat4 <- data.frame(w1, w2, L1_1, L2_1, a1, c, L1_2, L2_2, a2, y)
sum(dat4$c==1)/2000 #0.122

dat4$L1_2[dat4$c==1] <- NA
dat4$L2_2[dat4$c==1] <- NA
dat4$a2[dat4$c==1] <- NA
dat4$y[dat4$c==1] <- NA

dat4$c[dat4$c==0] <- "uncensored"
dat4$c[dat4$c==1] <- "censored"
dat4$c <- as.factor(dat4$c)

complete_dat4 <- dat4[complete.cases(dat4), ]

mytrue <- function(a1s = 1, a2s = 1){
  #visit 1;
  w1sim <- rbinom(samplesize, 1, prob = 0.5) #50-50 male and female;
  w2sim <- rnorm(samplesize, mean = 12, sd = 4) #age;
  L1_1sim <- rbinom(samplesize, 1, prob = plogis(-1 + 0.01*w1sim + 0.01*w2sim)) #a age related baseline binary clinical variable;
  L2_1sim <- rnorm(samplesize, mean = (0.01*w1sim + 0.01*w2sim), sd = 1) #a age & sex related baseline continuous clinical variable;
  
  #Visit 2, simulate all potential variables;
  L1_2sim <- rbinom(samplesize, 1, plogis(-1 + 0.01*w1sim + 0.01*w2sim + 0.2*L1_1sim - 0.2*a1s))
  L2_2sim <- rnorm(samplesize, mean = (L2_1sim - 0.01*w1sim + 0.01*w2sim - 0.2*a1s), sd = 1)
  
  #Visit 3, simulate potential outcomes;
  truemean = (1 - 1*a1s - 2*a2s + 0.01*a1s*a2s - 0.1*w1sim + 0.01*w2sim + 0.1*L1_2sim + 0.1*L2_2sim)
  return(mean(truemean))
}

mytrue(a1s = 1, a2s = 1) #potential outcome Y(a1=1, a2=1);
mytrue(a1s = 0, a2s = 0)
mytrue(a1s = 1, a2s = 0)
mytrue(a1s = 0, a2s = 1)
mytrue(a1s = 1, a2s = 1) - mytrue(a1s = 0, a2s = 0); #-3.010579, this can change due to random sampling;

library(arm)
library(ltmle)
library(ipw) #MSM;
library(gfoRmula) #g-computation;

# tmle without superlearner + kitchen sink gform + null Qform (kitchen sink);
tmle_model <- ltmle(dat4,
                    Anodes = c ("a1","a2"),
                    Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                    Ynodes = c("y"), 
                    Cnodes = c("c"),
                    survivalOutcome =FALSE,
                    # Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a1",
                    #            y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a1 + a2 + a1*a2"),
                    gform = c("a1 ~ w1 + w2 + L1_1 + L2_1",
                              "c ~ w1 + w2 + L1_1 + L2_1 +a1",
                              "a2 ~ w1 + w2 + L1_1 + L2_1 + a1 + L1_2 + L2_2"),
                    # gcomp = TRUE,
                    # iptw.only = FALSE,
                    # variance.method = "tmle",
                    # SL.library = "default", #with superlearner;
                    abar = list(c(1,1), c(0,0)),
                    estimate.time = FALSE)


summary(tmle_model, estimator="tmle")

# tmle with superlearner + kitchen sink gform + null Qform;
tmle_model <- ltmle(dat4,
                    Anodes = c ("a1","a2") ,
                    Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                    Ynodes = c("y"), 
                    Cnodes = c("c"),
                    survivalOutcome =FALSE,
                    # Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a1",
                    #            y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a1 + a2 + a1*a2"),
                    gform = c("a1 ~ w1 + w2 + L1_1 + L2_1",
                              "c ~ w1 + w2 + L1_1 + L2_1 +a1",
                              "a2 ~ w1 + w2 + L1_1 + L2_1 + a1 + L1_2 + L2_2"),
                    # gcomp = TRUE,
                    # iptw.only = FALSE,
                    # variance.method = "tmle",
                    SL.library = "default", #with superlearner;
                    abar = list(c(1,1), c(0,0)),
                    estimate.time = FALSE)

summary(tmle_model, estimator="tmle")

#tmle package for gcomputation + kitchen sink gform + null Qform;
tmle_model <- ltmle(dat4,
                    Anodes = c ("a1","a2") ,
                    Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                    Ynodes = c("y"), 
                    Cnodes = c("c"),
                    survivalOutcome =FALSE,
                    # Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a1",
                    #            y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a1 + a2 + a1*a2"),
                    gform = c("a1 ~ w1 + w2 + L1_1 + L2_1",
                              "c ~ w1 + w2 + L1_1 + L2_1 +a1",
                              "a2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a1"),
                    gcomp = TRUE,
                    # iptw.only = FALSE,
                    # variance.method = "tmle",
                    abar = list(c(1,1), c(0,0)),
                    estimate.time = FALSE)

summary(tmle_model)

#tmle package for gcomputation + kitchen sink gform + correct Qform;
tmle_model <- ltmle(dat4,
                    Anodes = c ("a1","a2") ,
                    Lnodes = c ("L1_1", "L2_1", "L1_2", "L2_2"), 
                    Ynodes = c("y"), 
                    Cnodes = c("c"),
                    survivalOutcome =FALSE,
                    Qform = c( L1_2 = "Q.kplus1 ~ w1 + w2 + L1_1 + a1",
                               y = "Q.kplus1 ~ w1 + w2 + L1_2 + L2_2 + a1 + a2 + a1*a2"),
                    gform = c("a1 ~ w1 + w2 + L1_1 + L2_1",
                              "c ~ w1 + w2 + L1_1 + L2_1 +a1",
                              "a2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a1"),
                    gcomp = TRUE,
                    # iptw.only = FALSE,
                    # variance.method = "tmle",
                    abar = list(c(1,1), c(0,0)),
                    estimate.time = FALSE)

summary(tmle_model)
