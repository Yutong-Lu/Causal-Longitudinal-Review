

data {
  int<lower=0> N; // Number of observations
  int<lower=0> N_cov; // Number of covariates
  matrix[N, N_cov] x1; // Covariates at visit 1
  matrix[N, N_cov] x2; // Covariates at visit 2
  int z1[N]; // treatment at visit 1 - conditional
  int z2[N]; // treatment at visit 2
  int z1_s[N]; // marginal
  int z2_s[N]; 
}

parameters {
  real alpha1; // Conditional treatment model intercept - visit 1
  real alpha2; // Conditional treatment model intercept - visit 2
  vector[N_cov] beta1; // Conditional treatment model covariate effects (5 covariates) - visit 1
  vector[N_cov] beta2; // Conditional treatment model covariate effects (5 covariates) - visit 2
  real beta_hist; // Conditional treatment model historical treatment effect (visit 1 treatment)
  
  real alpha1_s; // Marginal treatment model intercept - visit 1
  real alpha2_s; // Marginal treatment model intercept - visit 2
  real beta_hist_s; // Marginal treatment model historical treatment effect (visit 1 treatment)

  real<lower=0> sigma;
 
}


model {
  for (i in 1:N){
    z1[i] ~ bernoulli_logit(alpha1 + dot_product(beta1[1:N_cov], x1[i,]));
    z2[i] ~ bernoulli_logit(alpha2 + dot_product(beta2[1:N_cov], x2[i,])  + beta_hist * z1[i]);
    z1_s[i] ~ bernoulli_logit(alpha1_s);
    z2_s[i] ~ bernoulli_logit(alpha2_s + beta_hist_s * z1_s[i]); // try moving p's into model section of stan
  }
    

     // Priors
     alpha1 ~ normal(0, 5); // reducing variance of priors to 5
     alpha2 ~ normal(0, 5);
     beta1[1:N_cov] ~ normal(0, 5);
     beta2[1:N_cov] ~ normal(0, 5);
     beta_hist ~ normal(0, 5);
     
     alpha1_s ~ normal(0, 5);
     alpha2_s ~ normal(0, 5);
     beta_hist_s ~ normal(0, 5);
}

generated quantities{
  vector[N] weight_numerator;
  vector[N] weight_denominator;
  vector[N] z1_pred;
  vector[N] z2_pred;
  
  for(i in 1:N){
    weight_numerator[i] = alpha1_s*(alpha2_s + beta_hist_s * z1_s[i]);
    weight_denominator[i] = (alpha1 + dot_product(beta1[1:N_cov], x1[i,]))*(alpha2 + dot_product(beta2[1:N_cov], x2[i,])  + beta_hist * z1[i]);
     z1_pred[i] = bernoulli_logit_rng(alpha1 + dot_product(beta1[1:N_cov], x1[i,]));
     z2_pred[i] = bernoulli_logit_rng(alpha2 + dot_product(beta2[1:N_cov], x2[i,])  + beta_hist * z1[i]);
  }
    // take out weight_numerator (p1_s*p2_s). will save for each iteration, and then taking mean of all. 
    // take out weight_denominator (p1 * p2). Outside of this file, take mean of both weight numerator+denominator (colSums/colSums)
    // Ncolumns = # of patients, Nrows = # of iterations
}


