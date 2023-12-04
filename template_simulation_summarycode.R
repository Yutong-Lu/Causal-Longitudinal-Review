
setwd("C:/Users/kuan liu/Desktop/myrun/paper2/feb11_results")

library(xtable)

#1. h, 5, 300, mc;
results1<-read.table("paper2feb7_2005q1ec_1.txt")
results1<-results1[,-1]
#2. h, 5, 300, hc;
results2<-read.table("paper2feb7_2005q1hc_1.txt")
results2<-results2[,-1]
#3. l, 5, 300, mc;
results3<-read.table("paper2feb7_2005q3ec_1.txt")
results3<-results3[,-1]
#4. l, 5, 300, hc;
results4<-read.table("paper2feb7_2005q3hc_1.txt")
results4<-results4[,-1]



results<-rbind(results1,results2, results3, results4)
# results<-results[,-1]


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
tab <- rbind(tab,tablerow(estimator=c('MSM'), pevalues=results[,5], vevalues=results[,9], true=mean(results[,1]), cilvalues=xxx, ciuvalues=xxx))
tab <- rbind(tab,tablerow('g-comp', results[,6], results[,10], mean(results[,2])))
tab <- rbind(tab, tablerow('tlme', results[,7], results[,11], mean(results[,3])))
tab <- rbind(tab, tablerow('tlmesuper', results[,8], results[,12], mean(results[,4])))

for (i in 1:ncol(tab))
  if (is.numeric(tab[,i]))
    tab[,i] <- round(tab[,i], 3)
#tab[1,6]<-NA
print(tab)

xtable(tab)
