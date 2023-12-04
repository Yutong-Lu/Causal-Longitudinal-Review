# Weightit and survey package
# example code for frequentist MSM

# specify our weight models
Wmsm.out <- weightitMSM(list(z1 ~ L1_1 + L2_1,
                             z2 ~ L1_1 + L2_1 + L1_2 + L2_2 + Z1),
                        data = dat1, method = "ps",
                        stabilize = TRUE)

# take a look at the quality of the weights 
summary(Wmsm.out) #good to check for extreme weights (>10);


# estimate treatment effect
# create a survey object that uses the MSM weights 
# as survey weight to be modelled in the weighted GLM;
d.w.msm <- svydesign(~1, weights = Wmsm.out$weights, data = dat1)

#weighted GLM;
cum.fit <- svyglm(y ~ z1*z2, design = d.w.msm)
summary(cum.fit)