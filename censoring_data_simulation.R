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


# Visit 1
w1 <- rbinom(ntot, 1, prob = 0.5) #50-50 male and female (reference);
w2 <- rnorm(ntot, mean = 12, sd = 4) #age;
L1_1 <- rbinom(ntot, 1, prob = expit(-1 + 0.01*w1 + 0.01*w2)) #a age & sex related baseline binary clinical variable;
L2_1 <- rnorm(ntot, mean = (0.01*w1 + 0.01*w2), sd = 1) #a age & sex related baseline continuous clinical variable;

a_1 <- rbinom(ntot, 1, expit(-1 - 0.01*w1 + 0.01*w2 - 0.1*L1_1 + 0.1*L2_1))

# observational data;
L1_2 <- rbinom(ntot, 1, expit(-1 + 0.01*w1 + 0.01*w2 + 0.2*L1_1 - 0.2*a_1))
L2_2 <- rnorm(ntot, mean = (L2_1 - 0.01*w1 + 0.01*w2 - 0.2*a_1), sd = 1)
a_2 <- rbinom(ntot, 1, expit(-1 - 0.01*w1 + 0.01*w2 - 0.1*L1_2 + 0.1*L2_2 + a_1))

# end-of-study outcome;
y <- rnorm(ntot, mean = (1 - 1*a_1 - 2*a_2 + 0.01*a_1*a_2 - 0.2*w1 + 0.1*w2 + 0.4*L1_2 + 0.2*L2_2) , sd = 1)

# simulate c from a binary + logistic format;

c <- rbinom(ntot, 1,
              prob = expit(-3 + 0.002*w1 + 0.001*w2 + 0.003*L1_1 + 0.005*L2_1 + 
                               0.01*a_1))

dat4 <- data.frame(w1, w2, L1_1, L2_1, a_1, L1_2, L2_2, c, a_2, y)

dat4$y[dat4$c == 1] <- NA
# dat4$L1_2[dat4$c == 1] <- NA
# dat4$L2_2[dat4$c == 1] <- NA
dat4$a_2[dat4$c == 1] <- NA

sum(dat4$c)/1000

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
  truemean = (1 - 1*a_1s - 2*a_2s + 0.01*a_1s*a_2s - 0.1*w1sim + 0.01*w2sim + 0.1*L1_2sim + 0.1*L2_2sim)
  return(mean(truemean))
}

mytrue(a_1s = 1, a_2s = 1) #potential outcome Y(a_1=1, a_2=1);
mytrue(a_1s = 0, a_2s = 0)
mytrue(a_1s = 1, a_2s = 0)
mytrue(a_1s = 0, a_2s = 1)
mytrue(a_1s = 1, a_2s = 1) - mytrue(a_1s = 0, a_2s = 0); #-3.010579, this can change due to random sampling;

# Saving as CSV
write_csv(dat4,file = "datasets/censoring_data.csv")