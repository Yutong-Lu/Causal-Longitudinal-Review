###
#
###

library(tidyverse)

#pre-define function;
expit <- function(x){
  x <- exp(x)/(exp(x)+1) 
  return(x)
}

# Simulation setup;
set.seed(123)
ntot=1000;
samplesize=10000


# Visit 1;
w1 <- rbinom(ntot, 1, prob = 0.5) #50-50 male and female (reference);
w2 <- rnorm(ntot, mean = 12, sd = 4) #age;
L1_1 <- rbinom(ntot, 1, prob = expit(-1 + 0.01*w1 + 0.01*w2)) #a age & sex related baseline binary clinical variable;
L2_1 <- rnorm(ntot, mean = (0.01*w1 + 0.01*w2), sd = 1) #a age & sex related baseline continuous clinical variable;

a_1 <- rbinom(ntot, 1, expit(-1 - 0.01*w1 + 0.01*w2 - 0.1*L1_1 + 0.1*L2_1)); #exposure;

# observational data;
L1_2 <- rbinom(ntot, 1, expit(-1 + 0.01*w1 + 0.01*w2 + 0.2*L1_1 - 0.2*a_1))
L2_2 <- rnorm(ntot, mean = (L2_1 - 0.01*w1 + 0.01*w2 - 0.2*a_1), sd = 1)
a_2 <- rbinom(ntot, 1, expit(-1 - 0.01*w1 + 0.01*w2 - 0.1*L1_2 + 0.1*L2_2 + a_1))
# table(a_1, a_2)

# end-of-study outcome;
y <- rbinom(ntot, 1, prob = expit(- 1 + log(0.9)*a_1 + log(0.85)*a_2 + log(1.01)*a_1*a_2 
                                  + log(0.95)*w1 + log(1.05)*w2 + log(1.15)*L1_2 + log(1.1)*L2_2))

# look at simulating binary from odds ratios, different ways to do it
# table(y);
# sum(y)/ntot; #proportional of patients with a positive binary outcome (1);

# saving final data;
dat2 <- data.frame(w1, w2, L1_1, L2_1, a_1, L1_2, L2_2, a_2, y)

mytrue <- function(a_1s = 1, a_2s = 1){
  #visit 1;
  w1sim <- rbinom(samplesize, 1, prob = 0.5) #50-50 male and female;
  w2sim <- rnorm(samplesize, mean = 12, sd = 4) #age;
  L1_1sim <- rbinom(samplesize, 1, prob = expit(-1 + 0.01*w1sim + 0.01*w2sim)) #a age related baseline binary clinical variable;
  L2_1sim <- rnorm(samplesize, mean = (0.01*w1sim + 0.01*w2sim), sd = 1) #a age & sex related baseline continuous clinical variable;
  
  #Visit 2, simulate all potential variables;
  L1_2sim <- rbinom(samplesize, 1, expit(-1 + 0.01*w1sim + 0.01*w2sim + 0.2*L1_1sim - 0.2*a_1s))
  L2_2sim <- rnorm(samplesize, mean = (L2_1sim - 0.01*w1sim + 0.01*w2sim - 0.2*a_1s), sd = 1)
  
  #Visit 3, simulate potential outcomes;
  # calc odds of y=1
  true_prob = expit(- 1 + log(0.9)*a_1s + log(0.85)*a_2s + log(1.01)*a_1s*a_2s 
                   + log(0.95)*w1sim + log(1.05)*w2sim + log(1.15)*L1_2sim + log(1.1)*L2_2sim)
  return(mean(true_prob))
}

# convert prob into odds
# N0TE: plogis and expit are very close but not exactly the same
prob_11<-mytrue(a_1s = 1, a_2s = 1) #potential outcome Y(a_1=1, a_2=1);
prob_00<-mytrue(a_1s = 0, a_2s = 0)
prob_10<-mytrue(a_1s = 1, a_2s = 0)
prob_01<-mytrue(a_1s = 0, a_2s = 1)

# True absolute risk difference, -0.06414264; -0.06625605
mean(prob_11) - mean(prob_00)
# True risk ratio (relative risk), 0.8425255; 0.8383083
mean(prob_11)/mean(prob_00)
# True odds ratio, 0.7602477; 0.7537021
mean(prob_11/(1-prob_11))/mean(prob_00/(1-prob_00))


#alternative way to calculate true potential outcome;
# potential variables at visit 2 given treatment;
# L1_2_0 <- rbinom(ntot, 1, expit(-1 + 0.01*w1 + 0.01*w2 + 0.2*L1_1))
# L1_2_1 <- rbinom(ntot, 1, expit(-1 + 0.01*w1 + 0.01*w2 + 0.2*L1_1- 0.2))
# L2_2_0 <- rnorm(ntot, mean = (L2_1 - 0.01*w1 + 0.01*w2), sd = 1)
# L2_2_1 <- rnorm(ntot, mean = (L2_1 - 0.01*w1 + 0.01*w2 - 0.2), sd = 1)

# potential outcomes and the true simulated causal contrast;
# prob_00 <- expit(- 1 - 0.2*w1 + 0.1*w2 + 0.4*L1_2_0 + 0.2*L2_2_0)
# prob_10 <- expit(- 1 + log(0.9) - 0.2*w1 + 0.1*w2 + 0.4*L1_2_1 + 0.2*L2_2_1)
# prob_01 <- expit(- 1 + log(0.85) - 0.2*w1 + 0.1*w2 + 0.4*L1_2_0 + 0.2*L2_2_0)
# prob_11 <- expit(- 1 + log(0.9) + log(0.85) + log(1.01) - 0.2*w1 + 0.1*w2 + 0.4*L1_2_1 + 0.2*L2_2_1)

# # True absolute risk difference, -0.07372847;
# mean(prob_11) - mean(prob_00)
# # True risk ratio (relative risk), 0.8695818;
# mean(prob_11)/mean(prob_00)
# # True odds ratio, 0.7268491;
# mean(prob_11/(1-prob_11))/mean(prob_00/(1-prob_00))

# Saving as CSV
write_csv(dat2,file = "datasets/binary_outcome_data.csv")


