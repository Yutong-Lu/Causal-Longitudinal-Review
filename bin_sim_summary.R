
setwd("/Users/lorrainelu/Documents/GitHub/longitudinal_causal_analysis")

library(xtable)

results<-read.table("Nov30binrun.txt")
results<-results[,-c(1,2)]

results1<-read.table("Dec1binrun.txt")[,2]
results<-cbind(results1, results)

tablerow <- function(estimator, pevalues, vevalues, truevalue, cilvalues=NULL, ciuvalues=NULL) {
  return(data.frame(Estimator=estimator,
                    Mean=mean(pevalues),
                    RB=100.0 * mean((pevalues - truevalue)/abs(truevalue)),
                    # RB2=100.0 * median((pevalues - truevalue)/abs(truevalue)),
                    SD=sd(pevalues),
                    SE=mean(vevalues),
                    SERB=100.0 * (mean(vevalues) - sd(pevalues))/sd(pevalues),
                    CP=ifelse(is.null(cilvalues),
                              100.00 * mean((truevalue > (pevalues + qnorm(0.025) * vevalues)) & (truevalue < (pevalues + qnorm(0.975) * vevalues))),
                              100.00 * mean((truevalue > cilvalues) & (truevalue < ciuvalues)))
                    # ,CP_L=ifelse(is.null(cilvalues),
                    #             mean((pevalues) + qnorm(0.025) * sqrt(vevalues)),
                    #             mean(cilvalues)),
                    # CP_U=ifelse(is.null(cilvalues),
                    #             mean(pevalues) + qnorm(0.975) * mean(sqrt(vevalues)),
                    #             mean(ciuvalues))
  ))
}


tab <- NULL
tab <- rbind(tab,tablerow(estimator=c('MSM:simulation'), pevalues=results[,14], 
                          vevalues=results[,15], truevalue=results[,1], 
                          cilvalues=results[,16], ciuvalues=results[,17]))
tab <- rbind(tab,tablerow(estimator=c('MSM:kitchen sink'), pevalues=results[,18], 
                          vevalues=results[,19], truevalue=results[,1], 
                          cilvalues=results[,20], ciuvalues=results[,21]))
tab <- rbind(tab,tablerow('g-comp', results[,2], results[,3], results[,1], 
                          cilvalues=results[,4], ciuvalues=results[,5]))
tab <- rbind(tab, tablerow('tlme', results[,6], results[,7], results[,1],
                           cilvalues=results[,8], ciuvalues=results[,9]))
tab <- rbind(tab, tablerow('tlmesuper', results[,10], results[,11], results[,1],
                           cilvalues=results[,12], ciuvalues=results[,13]))

for (i in 1:ncol(tab))
  if (is.numeric(tab[,i]))
    tab[,i] <- round(tab[,i], 5)
#tab[1,6]<-NA
print(tab)
# 
# xtable(tab)

mode <- function(x){
  uniqv <- unique(x)
  uniqv[which.max(tabulate(match(x, uniqv)))]
}
mode(results[,1])

plot(density(results[,1]))
abline(v=0.78)
abline(v=0.765)
