### JAGS model
cat( " model {

     #N = nobs
     for (i in 1:N) {

     # conditional treatment assignment model;

     z1[i] ~ dbern(p1[i])
     logit(p1[i]) <- b10 + x1[i,] %*% b11

     z2[i] ~ dbern(p2[i])
     logit(p2[i]) <- b20 + x2[i,] %*% b21 + b_hist*z1[i]

     # marginal treatment assignment model;

     z1s[i] ~ dbern(p1s[i])
     logit(p1s[i]) <- bs10

     z2s[i] ~ dbern(p2s[i])
     logit(p2s[i]) <- bs20 + bs_hist*z1[i]
     }

     # Priors
     b10 ~ dnorm(0,5) #true 0;
     
     for(j in 1:Ncov){
        b11[j] ~ dnorm(0,5) #true 0.2;
        b21[j] ~ dnorm(0,5)
     }
     
     b20 ~ dnorm(0,5)
     b_hist ~ dnorm(0,5)

     bs10 ~ dnorm(0,5) #true 0;
     bs20 ~ dnorm(0,5)
     bs_hist ~ dnorm(0,5)

     }",
     file = "model_bmsm.txt")
