library(SimDesign)
library(rstan)
library(tidyverse)
library(gfoRmula)
library(ltmle)
library(gtools)
library(WeightIt)
library(survey)
library(R2jags)
library(runjags)
library(bayesplot)
library(rstanarm)

expit<-function(x){1/(1+exp(-x))}
logit<-function(p){log(p)-log(1-p)}


loglik2<-function(param,resp,mmat, weight){
  
  n<-length(resp)
  ### NOTE: param will have 5 entries, one for each theta (intercept, visit 1 treatment effect, visit 2 treatment effect, visit 1-2 interaction effect) and sigma (response)
  theta <- param[1:4] 
  sigma <- param[5]  ## only 1 sigma because end of study outcome.
  e<- (resp - mmat%*%theta)
  logl <- NA
  
  if (sigma > 0){
    logl1<- -0.5*log(sigma^2) - ((e)^2)/(2*(sigma^2)) 
    logl<- sum(weight*logl1)
  }
  return(logl)
  
}


mysim <- function(from=1, to=1, ntot=500, ncov=5, ncov_d = 2,
                  mean_x1=-0.8, mean_x2=0.5, mean_x3=-0.8, mean_x4=0.8, mean_x5=0.5, mean_y=0.8, 
                  sigma_x1 = 1, sigma_x2 = 1, sigma_x3=1, sigma_x4=1, sigma_x5=1, sigma_y = 1,
                  rho12 = -0.5, rho13 = 0.5, rho14 = 0.5, rho15=-0.5,
                  rho23 = -0.5, rho24 = -0.5, rho25 = 0.5,
                  rho34 = 0.5, rho35 = -0.5, rho45 = -0.5,
                  coeffz_x1 = 0.1, coeffz_x2 = -0.1, coeffz_x3 = 0.02, coeffz_x4 = 0.02, coeffz_x5 = -0.1, coeffz_z = 1,
                  coeffx1_z = 0.1, coeffx2_z = -0.1, coeffx3_z = 0.01, coeffx4_z = 0.01, coeffx5_z = -0.1,
                  coeffy_x1 = 0.01, coeffy_x2 = -0.01, coeffy_x3 = 0.01, coeffy_x4 = 0.01, coeffy_x5 = -0.01,
                  nobs = 500,
                  alpha = 5,
                  effect = 0.5,
                  lag_effect = 0.5^2,
                  effect_int=-0.3,
                  niter=15000, nburnin=5000, thin=5, nboot=2000) {
  # 
  # i = 1;
  # # Range of Number of simulation runs - lower bound
  # from=1;
  # # Range of Number of simulation runs - upper bound
  # to=1;
  # # Sample size
  # ntot=500;
  # # Number of covariates
  # ncov=5;
  # 
  # # Number of discrete variables
  # ncov_d=2;
  # 
  # # Mean of X's and Y
  # mean_x1=-0.8; mean_x2=0.5; mean_x3=-0.8; mean_x4=0.8; mean_x5=0.5; mean_y=0.8;
  # # Standard deviances of X's and Y
  # sigma_x1 = 1; sigma_x2 = 1; sigma_x3=1; sigma_x4=1; sigma_x5=1; sigma_y = 1;
  # 
  # # Correlation coefficients between pairwise X's
  # rho12 = -0.5; rho13 = 0.5; rho14 = 0.5; rho15=-0.5;
  # rho23 = -0.5; rho24 = -0.5; rho25 = 0.5;
  # rho34 = 0.5; rho35 = -0.5; rho45 = -0.5;
  # 
  # # Treatment assignment model coefficients
  # coeffz_x1 = 0.1; coeffz_x2 = -0.1; coeffz_x3 = 0.02; coeffz_x4 = 0.02; coeffz_x5 = -0.1; coeffz_z = 1;
  # 
  # # Treatment effects on X's
  # coeffx1_z = 0.1; coeffx2_z = -0.1; coeffx3_z = 0.01; coeffx4_z = 0.01; coeffx5_z = -0.1;
  # 
  # # Covariate effects on Y
  # coeffy_x1 = 0.01; coeffy_x2 = -0.01; coeffy_x3 = 0.01; coeffy_x4 = 0.01; coeffy_x5 = -0.01;
  # # coeffy_x = matrix(c(coeffy_x1, coeffy_x2, coeffy_x3, coeffy_x4, coeffy_x5), ncov, 1);
  # 
  # # treatment effects on Y;
  # alpha = 5;
  # effect = 0.5;
  # lag_effect = 0.5^2;
  # effect_int=-0.3;
  # # 
  # Empty y
  y =array(NA, c(2, 2, ntot));


  coeffz_x = matrix(c(coeffz_x1, coeffz_x2, coeffz_x3, coeffz_x4, coeffz_x5), ncov, 1);
  coeffy_x = matrix(c(coeffy_x1, coeffy_x2, coeffy_x3, coeffy_x4, coeffy_x5), ncov, 1);
  niter=15000; nburnin=5000; thin=5; nboot=2000;

    samplesize <- (niter - nburnin)/thin
  results <- matrix(NA, to, 26)
  results_names<- matrix(NA, to, 26)
  
  ## Number of variables to keep track of
  nL <-ncov
  
  
  for (i in from:to) {
    
    ## Setting different seed for each simulation run
    iter <- 66500 + (i - 1)
    set.seed(iter)
    
    ## Matrix tracking all covariate values for all 2 visits, under each of the 2 Potential Outcomes for the first 2 visits. 
    #first dimension is 2; represents z1 = 0 and z1 = 1; 
    #second dimension is 2; represents z2 = 0 and z2 = 1; 
    #third dimension is ntot; represents number of subjects; 
    #fourth dimension k*nL represent visits for x1-x5, and x1d and x2d (the discrete versions of x1 and x2), in this order;
    L <- array(NA, c(2, 2, ntot, 2 * (nL))) 
    
    ## Matrix storing "observed" covariate values, so there is only 1 Potential Outcome per visit.
    Lobs <- matrix(NA, ntot, 2 * (nL))
    
    ## Observed treatment assignment
    z <- matrix(NA, ntot, 2)
    
    
    ### 1. visit 1:
    
    ## cov matrix of x1,x2,x3,x4,x5 
    s <- matrix(c(sigma_x1^2, (rho12)*sigma_x1*sigma_x2, (rho13)*sigma_x1*sigma_x3, (rho14)*sigma_x4*sigma_x1, (rho15)*sigma_x5*sigma_x1,
                  (rho12)*sigma_x1*sigma_x2, sigma_x2^2, (rho23)*sigma_x3*sigma_x2, (rho24)*sigma_x4*sigma_x2, (rho25)*sigma_x5*sigma_x2,
                  (rho13)*sigma_x1*sigma_x3, (rho23)*sigma_x2*sigma_x3, sigma_x3^2, (rho34)*sigma_x4*sigma_x3, (rho35)*sigma_x5*sigma_x3,
                  (rho14)*sigma_x1*sigma_x4, (rho24)*sigma_x2*sigma_x4, (rho34)*sigma_x3*sigma_x4, sigma_x4^2, (rho45)*sigma_x5*sigma_x4,
                  (rho15)*sigma_x1*sigma_x5, (rho25)*sigma_x2*sigma_x5, (rho35)*sigma_x3*sigma_x5, (rho45)*sigma_x4*sigma_x5, sigma_x5^2), 5, 5) 
    #test; 
    # det(s)
    
    ## Generating continuous X's
    Lobs[,1:nL] <- L[1,1,,1:nL] <- L[1,2,,1:(nL)] <- L[2,1,,1:(nL)] <- L[2,2,,1:(nL)] <- rmvnorm(ntot, c(mean_x1, mean_x2, mean_x3, mean_x4, mean_x5), s)
    ## Linear predictor of treatment assignment model
    # m_z<- Lobs[,1] * coeffz_x1 + Lobs[,2] * coeffz_x2 + Lobs[,3] * coeffz_x3 + Lobs[,4] * coeffz_x4 + Lobs[,5] * coeffz_x5 
    m_z<- Lobs[,1:5] %*% coeffz_x
    ## Treatment assignment probability
    pz <- expit(m_z)
    ## Observed treatment assignment
    z[,1] <- rbinom(ntot,1,prob=pz)
    
    
    ### 2. visit 2:
    # for Z1 = 0;
    
    ## Generating continuous X's
    L[1,1,,(nL+1):(2*nL)] <- L[1,2,,(nL+1):(2*nL)] <- t(apply(L[1,1,,1:(nL)], 1, function(x){rmvnorm(1, mean = x, s)}))
    
    # for z1 = 1;
    
    ## Treatment effects on covariates
    trt<- matrix(c(rep(coeffx1_z,ntot), rep(coeffx2_z,ntot), rep(coeffx3_z,ntot), rep(coeffx4_z,ntot), rep(coeffx5_z,ntot)),ntot,ncov)
    
    ## Generating continuous X's
    L[2,1,,(nL+1):(2*nL)] <- L[2,2,,(nL+1):(2*nL)] <- t(apply( (L[2,1,,1:(nL)]) + trt, 1, function(x){rmvnorm(1, mean = x, s)})) 
    
    ## Extracting only the observed covariate values (based on observed treatment assignment)
    for (j in (nL+1):(2*nL)) {
      Lobs[,j] <- L[cbind(z[,1]+1,1,1:ntot,j)]
    }
    
    ## Linear predictor of treatment assignment model
    m_z<- Lobs[,(nL+1):(2*nL)] %*% coeffz_x + coeffz_z*z[,1]
    ## Treatment assignment probability
    pz <- expit(m_z)
    ## Observed treatment assignment
    z[,2] <- rbinom(ntot,1,prob=pz)
    
    
    ## Cumulative treatments
    zobs <- rowSums(z)
    table(zobs)
    
    ## Observe outcomes
    lp <- alpha + effect*z[,1] + (lag_effect)*z[,2] + effect_int*z[,1]*z[,2] + Lobs[,(nL + 1):(2*nL)]  %*% coeffy_x
    yobs <- rnorm(n = nrow(Lobs), mean = lp, sd = sigma_y)
    
    
    # Kuan : old code to get the true y and true causal effect; Outcomes
    ylong <- NULL
    zlong <- NULL
    x2long <- NULL
    counter <- 1
    for (j in 0:1) {
      for (k in 0:1) {
        # j = 0
        # k = 0
        ## Response model linear predictor
        lp <- alpha + effect*j + (effect**2)*k + effect_int*j*k + cbind(L[j+1, k+1,,(nL + 1):(2*nL)]) %*% coeffy_x    ## Should these covariates be used?
        if(any(is.na(lp))){
          L_na <- L[j+1,k+1,,(nL+1):(2*nL)]
          lp_na <- lp
          print(j)
          print(k)
        }
        ## Applying link function (identity link) to get response values
        y[j+1,k+1,] <- lp
        
        ## Changing response values to long format
        ylong <- c(ylong, y[cbind(j+1,k+1,1:ntot)])
        
        ## Changing treatment values to long format
        zlong <- c(zlong, rep(j+k, ntot))
        
        ## Changing covariates at last visit to long format
        x2long <- x2long %>% bind_rows(L[j+1, k+1,,(nL+1):(2*nL)] %>% as.data.frame() %>% mutate(z1 = j, z2 = k))
      }
    }
    
    ## Extracting observed response values
    summary(ylong)
    summary(yobs)
    
    ## Potential treatment values for visits 1 and 2, respectively, in long format
    z1long <- c(rep(0, ntot*2), rep(1, ntot*2))
    z2long <- c(rep(0, ntot), rep(1, ntot), rep(0, ntot), rep(1, ntot))
    
    ## True Model
    truemodel <- lm(ylong ~ z1long*z2long + as.matrix(x2long[,1:5])) #in analysis, can investigate an incorrect likelihood where z1 and z2 have same effect
    summary(truemodel)

    # truemodel2 <- lm(ylong ~ z1long*z2long) #in analysis, can investigate an incorrect likelihood where z1 and z2 have same effect
    # summary(truemodel2)
    # mean(predict(truemodel2, newdata = data.frame(z1long = rep(1,1000), z2long=rep(1,1000)))) - mean(predict(truemodel2, newdata = data.frame(z1long = rep(0,1000), z2long=rep(0,1000))))
    # 0.4981684
    
    ## True model results
    # mean of Y_00 and Y_11;
    y_00<-mean(predict(truemodel, newdata = cbind(z1long = 0, z2long=0, (x2long[,1:5]))))
    y_11<-mean(predict(truemodel, newdata = cbind(z1long = 1, z2long=1, (x2long[,1:5]))))
    results[i,1] <- y_11  
    results_names[i, 1] <- "truemodel11"
    results[i,2] <- y_00  
    results_names[i, 2] <- "truemodel00"
    results[i,3] <- y_11 - y_00
    results_names[i, 3] <- "truemodelATE"
    
 
    
    #creating observed wide datasets;
    obs<-data.frame(Lobs,z, z[,1], rowSums(z[,1:2]),rep(1,ntot),rep(2,ntot), yobs)
    colnames(obs) <- c(paste0("x1", 1:ncov), paste0("x2", 1:ncov),"z1", "z2","cumz1","cumz2", "visit1","visit2", "y")
    obs$id<-seq(1, ntot, by=1)
    
    ## treatment Models
    tmodel1 <- glm(z1 ~ x11+x12+x13+x14+x15, data = obs, family = binomial(link =logit)) #in analysis, can investigate an incorrect likelihood where z1 and z2 have same effect
    tmodel2 <- glm(z2 ~ x21+x22+x23+x24+x25+z1, data = obs, family = binomial(link=logit)) #in analysis, can investigate an incorrect likelihood where z1 and z2 have same effect
    smodel1 <- glm(z1 ~ 1, data = obs, family = binomial(link=logit)) #in analysis, can investigate an incorrect likelihood where z1 and z2 have same effect
    smodel2 <- glm(z2 ~ z1, data = obs, family = binomial(link=logit)) #in analysis, can investigate an incorrect likelihood where z1 and z2 have same effect
    
    
    #Creating observed long datasets to model weighted GEE (Kuan: update the creation of long data to incorporate the new observed y);
    lobs<-rbind( cbind(seq(1, ntot, by=1),Lobs[,1:nL],z[,1]),cbind(seq(1, ntot, by=1), Lobs[,(nL+1):(2*nL)],z[,2]))
    time<-c(rep(1, ntot), rep(2, ntot))
    lcumz<-c(z[,1], rowSums(z[,1:2]))
    obslong<-data.frame(lobs, time, lcumz)
    colnames(obslong)<-c("id", "x1", "x2", "x3", "x4", "x5", "z", "time", "lcumz")
    yobs_df <- data.frame(id = 1:ntot, y = yobs)
    obslong <- obslong %>% left_join(yobs_df, by = "id")
    obslong$y <- with(obslong, ifelse(time < 2, NA, y))
    # dim(obslong) #1500 rows, 6 cols
    obslong <- obslong[order(obslong$id),] 

    
    ## Confounded Model
    confmodel <- lm(y ~ z1*z2 + x21 + x22 + x23 + x24 + x25, data = obs)
    summary(confmodel)
    ## Confounded model results
    # mean of Y_00 and Y_11;
    y_conf_00<-mean(predict(confmodel, newdata = cbind(z1 = 0, z2=0, (obs[,6:10]))))
    y_conf_11<-mean(predict(confmodel, newdata = cbind(z1 = 1, z2=1, (obs[,6:10]))))
    results[i,4] <- y_conf_11  
    results_names[i, 4] <- "confmodel11"
    results[i,5] <- y_conf_00  
    results_names[i, 5] <- "confmodel00"
    results[i,6] <- y_conf_11 - y_conf_00
    results_names[i, 6] <- "confmodelATE"
    
    ## G-COMPUTATION
    print("G-Computation")
    # gformula;
    # this package is picky do not include variables that are not used in the regression;
    obslong_new <- obslong %>% 
                      select(id, x1, x2, x3, x4, x5, time, z, y) %>%
                      mutate(time = time - 1)

    id <- 'id'
    time_name <- 'time'
    # create one to for the package
    covnames <- c("x1", "x2", "x3", "x4", "x5", "z")
    outcome_name <- 'y'
    covtypes <- c('normal', 'normal', 'normal', 'normal', 'normal', 'binary')
    histories <- c(lagged) ## Not sure what this is doing
    histvars <- list(c("x1", "x2", "x3", "x4", "x5",'z'))
    
    covparams <- list(covmodels = c(x1 ~ lag1_x1 + + lag1_z ,
                                    x2 ~ x1 + lag1_x2 + lag1_z ,
                                    x3 ~  x1 + x2 + lag1_x3 + lag1_z ,
                                    x4 ~ x1 + x2 + x3 + lag1_x4 + lag1_z ,
                                    x5 ~ x1 + x2 + x3 + x4 + lag1_x5 + lag1_z ,
                                    z ~ x1 + x2  + x3  + x4  + x5 + lag1_z))

    ymodel <- y ~ lag1_z + z + lag1_z*z + x1 + x2  + x3  + x4  + x5
    
    intvars <- list('z', 'z')
    interventions <- list(list(c(static, rep(0, 2))),
                          list(c(static, rep(1, 2))))
    int_descript <- c('Never treat', 'Always treat')
    
    gform_cont_eof2 <- gformula_continuous_eof(
      obs_data = obslong_new,
      id = id,
      time_name = time_name,
      covnames =covnames,
      outcome_name = outcome_name, 
      covtypes = covtypes,
      covparams = covparams,  ymodel = ymodel,
      intvars = intvars, interventions = interventions,
      int_descript = int_descript, 
      ref_int = 1,
      histories = c(lagged), histvars = histvars,
      nsimul = 1000,
      nsamples = 1000, parallel = TRUE, ncores = 4,
      seed = 123)
    
    gform_cont_eof2
    ## used to be overestimating - now only slightly 0.4925430
    results[i, 7] <- gform_cont_eof2$result$`Mean difference`[3]
    results_names[i, 7] <- "G-Comp ATE"
    results[i, 8] <- gform_cont_eof2$result$`MD lower 95% CI`[3]
    results[i, 9] <- gform_cont_eof2$result$`MD upper 95% CI`[3]
    results_names[i, 8] <- "G-Comp ATE lower 95% CI"
    results_names[i, 9] <- "G-Comp ATE upper 95% CI"
    
    ## TMLE
    print("TMLE")
    # tmle usign ltmle packages; 
    obs <- obs %>% dplyr::select(x11, x12, x13, x14, x15, x21, x22, x23, x24, x25, z1, z2, y)
    # ltmle with superlearner + kitchen sink gform + Qform;
    tmle_model_s <- ltmle(obs,
                          Anodes = c ("z1","z2") ,
                          Lnodes = c ("x11", "x12","x13", "x14", "x15", "x21", "x22", "x23", "x24", "x25"), 
                          Ynodes = c("y"), 
                          survivalOutcome =FALSE,
                          Qform = c(
                            # x21 = "Q.kplus1 ~ x11 + x12 + x13 + x14 + x15 + z1",
                            y = "Q.kplus1 ~ x21 + x22 + x23 + x24 + x25 + z1 + z2 + z1*z2"),
                          gform = c("z1 ~ x11 + x12 + x13 + x14 + x15",
                                    "z2 ~ x21 + x22 + x23 + x24 + x25 + z1"),
                          SL.library = "default", #with superlearner;
                          abar = list(c(1,1), c(0,0)),
                          estimate.time = FALSE)
    
    tmle_summary <- summary(tmle_model_s, estimator="tmle")
    results[i, 10] <- tmle_summary$effect.measures$ATE$estimate
    results[i, 11] <- tmle_summary$effect.measures$ATE$CI[1]
    results[i, 12] <- tmle_summary$effect.measures$ATE$CI[2]
    results_names[i, 10] <- "TMLE ATE"
    results_names[i, 11] <- "TMLE ATE lower 95% CI"
    results_names[i, 12] <- "TMLE ATE upper 95% CI"
    
    ## FREQUENTIST MSM
    print("Frequentist MSM")
    # Start bootstrap here;
    boot.est <- rep(NA, 1000)
    for (j in 1:1000){
      boot.idx <- sample(1:ntot,size = ntot, replace = T)
      boot.data <- obs[boot.idx,]
      
      Wmsm.out <- weightitMSM(list(z1 ~ x11 + x12 + x13 + x14 + x15,
                                 z2 ~ x21 + x22 + x23 + x24 + x25+ z1),
                            data = boot.data, method = "ps",
                            stabilize = TRUE)
    
    # estimate treatment effect
    d.w.msm <- svydesign(~1, weights = Wmsm.out$weights,
                         data = boot.data)
    
    cum.fit <- svyglm(y ~ z1 + z2 + z1*z2 + x21 + x22 + x23 + x24 + x25, design = d.w.msm)
    boot.est[j]<- sum(cum.fit$coefficients[c("z1", "z2", "z1:z2")])
    }
    
    
    sd(boot.est)
    c(quantile(boot.est, prob = 0.025), quantile(boot.est, prob = 0.975)) 
    results[i, 13] <- mean(boot.est)
    results[i, 14] <- quantile(boot.est, prob = 0.025)
    results[i, 15] <- quantile(boot.est, prob = 0.975)
    results_names[i, 13] <- "Frequentist MSM ATE (bootstrap)"
    results_names[i, 14] <- "Frequentist MSM ATE lower 95% CI (bootstrap)"
    results_names[i, 15] <- "Frequentist MSM ATE upper 95% CI (bootstrap)"
    
    ## ^ BOOTSTRAP FOR CI's
    
    

    ## Bayesian MSM
    print("Bayesian MSM")
    ### Stan Modelling
    stan_data <- list(N = nrow(Lobs),
                      N_cov = ncov,
                      x1 = Lobs[,1:nL],
                      x2 = Lobs[,(nL+1):(2*nL)],
                      z1 = z[,1],
                      z2 = z[,2],
                      z1_s = z[,1],
                      z2_s = z[,2])
    
    mod <- stan(data = stan_data,
                file = ("JAGS model in stan_previous.stan"),
                iter = 6000,
                seed = 243)
    # summary(mod)[["summary"]]
    saveRDS(mod, "BMSM_mod_stan_previous.RDS")
    # mod <- readRDS("BMSM_mod_stan_previous.RDS")
    print("done")
    
    # mod2 <- stan(data = stan_data,
    #             file = ("JAGS model in stan.stan"),
    #             iter = 6000,
    #             seed = 243)
    # summary(mod2)[["summary"]]
    # saveRDS(mod2, "BMSM_mod_stan.RDS")
    # mod2 <- readRDS("BMSM_mod_stan.RDS")
    

    # mod_summary <- summary(mod2)
    # pairs(mod2, pars = c("alpha1", "beta1"))
    
    ## Extracting the weights
    mod_weight_numerator <- rstan::extract(mod)[[c("weight_numerator")]] %>% as.data.frame()
    mod_weight_denominator <- rstan::extract(mod)[[c("weight_denominator")]] %>% as.data.frame()
    weights_ratio <- colSums(mod_weight_numerator)/colSums(mod_weight_denominator)
    wmean_s <- c(weights_ratio)
    
    ## Bootstrap estimates
    bootest3<-numeric(nboot)

    for (j8 in 1:nboot) {
      alpha <- as.numeric(rdirichlet(1, rep(1.0, ntot)))
      inits2 <- c(1, 1, 1, 1, 1)
      ## y_end ~ normal distribution with mean and sigma_squared. Mean is modelled by theta0+theta1*z1+theta2*z2 (if there's interaction, then add theta3*z1*z2). sigma has to be > 0. 0 for theta's and 1 for sigma as initial points. need design matrix
      maxim <- optim(inits2, fn=loglik2,resp=yobs, mmat=cbind(1, z[,1], z[,2], z[,1] * z[,2]), 
                     weight=alpha*wmean_s,
                     # control=list(fnscale=-1), method='BFGS', hessian=F)
                     control=list(fnscale=-1), method = "Nelder-Mead", hessian=F)
      bootest3[j8] <- sum(maxim$par[2:4])
      
      
      if (j8 %% 50 == 0) {
        print(j8)
      }
    }
      results[i, 16] <- mean(bootest3)
      results[i, 17] <- quantile(bootest3, prob = 0.025)
      results[i, 18] <- quantile(bootest3, prob = 0.975)
      results_names[i, 16] <- "Bayesian MSM ATE (bootstrap)"
      results_names[i, 17] <- "Bayesian MSM ATE lower 95% CI (bootstrap)"
      results_names[i, 18] <- "Bayesian MSM ATE upper 95% CI (bootstrap)"
      colnames(results) <- results_names
      
    # }
    
    
    
    
  }
  saveRDS(results, "results_BMSM.RDS")
  results
}

results <- mysim()  
colnames(results[1:26]) <- results[27:52]




