An Evaluation of Confidence Intervals for a Cumulative Proportion to
Enable Decisions at Interim Reviews of Single Arm Clinical Trials
================
Isabelle R. Weir and Linda J. Harrison
Janyary 9, 2024

### load R packages:

``` r
library(ggplot2)
library(survival)
library(km.ci)
library(binom) 
library(bpcp)
```

### Simulation Function (simulateCIs)

The function takes inputs for sample size (n), fixed time point (tau),
maximum censoring time (cens_max), target cumulative proportion
(target_CP), alpha type 1 error rate (alpha), and number of replicates
for the simulation (nrep). 
The function uses a uniform censoring pattern from 0 to week 8 (see C1 defined in the code below). To replicate the linear upwards and linear downwards censoring patterns, we replace C1 with the following:\
Linear Upwards: C1 <- sqrt(64\*runif(n=0.8\*scenarios$n[i]))\
Linear Downwards: C1 <- 8-sqrt(64\*runif(n=0.8\*scenarios$n[i]))\
A user can toggle the code to use these patterns in the code below.

``` r
simulateCIs <- function(n, tau, cens_max, target_CP, alpha, nrep){

     method          <- c("Greenwood clog-log",
                          "Thomas-Grunkemeier",
                          "BPCP",
                          "BPCP Mid P",
                          "Rothman Wilson",
                          "Clopper-Pearson",
                          "Wilson")
     
     ## compile scenarios based on function inputs:
     scenarios        <- expand.grid(n=n, 
                                     tau=tau,
                                     alpha=alpha,
                                     #exp_rate=exp_rate,
                                     target_CP=target_CP,
                                     cens_max=cens_max,
                                     nrep=nrep)
     
     scenarios$exp_rate <- c(log(1-scenarios$target_CP)/(-scenarios$tau))
     scenarios$scenID   <- 1:nrow(scenarios)
     scenarios
     
     ## prepare empty results dataframe
     allresults <- data.frame(scenID=numeric(), 
                              ENt=numeric(),
                              EYt=numeric(),
                              n=numeric(),
                              alpha=numeric(),
                              CPt1=numeric(),
                              method=character(),
                              lberr=numeric(),
                              uberr=numeric(), stringsAsFactors = FALSE)
     
     for(i in scenarios$scenID){
          
          print(paste0("Scenario ", i, " of ", nrow(scenarios), "; n=", scenarios$n[i], "; targetCP=", scenarios$target_CP[i]))
          print(Sys.time())
          
          result_i <- data.frame(scenID=numeric(length(method)), 
                                 n=numeric(length(method)),
                                 alpha=numeric(length(method)),
                                 CPt1=numeric(length(method)),
                                 cProp=numeric(length(method)),
                                 sProp=numeric(length(method)),
                                 method=character(length(method)),
                                 lberr=numeric(length(method)),
                                 uberr=numeric(length(method)), stringsAsFactors = FALSE)
          
          columns <- c("cProp", "sProp",
                       "lb_GW", "ub_GW", 
                       "lb_TG", "ub_TG",
                       "lb_BPCP", "ub_BPCP", 
                       "lb_BPCPmidp", "ub_BPCPmidp", 
                       "lb_RW", "ub_RW", 
                       "lb_Clopper", "ub_Clopper", 
                       "lb_W", "ub_W") 
          
          temp <- data.frame(matrix(nrow=nrep, ncol=length(columns)))
          colnames(temp) <- columns
          
          St1 <- exp(-(scenarios$exp_rate[i])*scenarios$tau[i])
     
          for(r in 1:nrep){
               
               # data generation with censoring
               E       <- rexp(n=scenarios$n[i], rate=scenarios$exp_rate[i])
               # Uniform Censoring
               C1      <- runif(n=0.8*scenarios$n[i], 0, 8)
               # Linear Upwards Censoring 
               # C1    <- sqrt(64*runif(n=0.8*scenarios$n[i]))
               # Linear Downwards Censoring 
               # C1    <- 8-sqrt(64*runif(n=0.8*scenarios$n[i]))

               C2      <- runif(n=0.2*scenarios$n[i], 8, 10)
               C       <- c(C1, C2)
               time   <- pmin(E, C)
               status <- as.numeric(E<=C)
               simdat <- data.frame(time=time, event=status, followup=C)
               
               ### CONFIDENCE INTERVALS FOR CUMULATIVE PROBABILITY
               
               # Greenwood clog-log 
               ## if no events before tau, then CI for cum prob is (0,0)
               ## if all events occur before tau, then CI for cum prob is (1,1)
               KM                  <- survfit(Surv(time, event)~1, data=simdat, 
                                              conf.int=1-scenarios$alpha[i], 
                                              conf.type="log-log")
               
               temp$lb_GW[r]       <- ifelse(sum(simdat[which(simdat$time<scenarios$tau[i]),]$event)==0, 0,
                                             ifelse(max(KM$time) < scenarios$tau[i], 1, 
                                                    1-summary(KM, times=scenarios$tau[i])$upper ))
               
               temp$ub_GW[r]       <- ifelse(sum(simdat[which(simdat$time<scenarios$tau[i]),]$event)==0, 0,
                                             ifelse(max(KM$time) < scenarios$tau[i], 1, 
                                                    1-summary(KM, times=scenarios$tau[i])$lower ))
               
               # save the cumulative proportion estimate
               temp$cProp[r] <- ifelse( max(KM$time) < scenarios$tau[i], 1-min(summary(KM)$surv), 
                                        1-summary(KM, times=scenarios$tau[i])$surv)
               
               ## Rothman-Wilson Score method
               ## if no events, then CI for cum prob is (0,0) automatically
               ## if ALL events occur before tau, then CI for cum prob is (1,1)
               KM2                 <- survfit(Surv(time, event)~1, conf.int=1-scenarios$alpha[i], data=simdat)
               
               RW                  <- km.ci(KM2, conf.level=1-scenarios$alpha[i], method="rothman")
               
               temp$lb_RW[r]       <- ifelse(max(RW$time) < scenarios$tau[i], 1, 
                                             1-summary(RW, times=scenarios$tau[i])$upper )
               
               temp$ub_RW[r]       <- ifelse(max(RW$time) < scenarios$tau[i], 1, 
                                             1-summary(RW, times=scenarios$tau[i])$lower )
               
               ## Thomas-Grunkemeier
               ## if no events, then CI for cum prob is (0,0) automatically with kmci1TG
               ## if ALL events occur before tau, then CI for cum prob is (1,1)
               TG                 <- kmci1TG(simdat$time, simdat$event, alpha=scenarios$alpha[i], tstar=scenarios$tau[i])
               
               temp$lb_TG[r]      <- ifelse(max(TG$time) < scenarios$tau[i], 1, 
                                            1-summary(TG, times=scenarios$tau[i])$upper)
               
               temp$ub_TG[r]      <- ifelse(max(TG$time) < scenarios$tau[i], 1, 
                                            1-summary(TG, times=scenarios$tau[i])$lower)
               
               ## BPCP (method of moments)
               BPCP                <- bpcpfit(simdat$time, simdat$event, conf.level=1-scenarios$alpha[i], nmc=0, midp=FALSE)
               temp$lb_BPCP[r]     <- 1-BPCP[[1]]$upper[which.max(BPCP[[1]]$R>scenarios$tau[i])]
               temp$ub_BPCP[r]     <- 1-BPCP[[1]]$lower[which.max(BPCP[[1]]$R>scenarios$tau[i])]
               
               ## BPCP Mid P
               BPCPmidp            <- bpcpfit(simdat$time, simdat$event, conf.level=1-scenarios$alpha[i], nmc=0, midp=TRUE)          
               temp$lb_BPCPmidp[r] <- 1-BPCPmidp[[1]]$upper[which.max(BPCPmidp[[1]]$R>scenarios$tau[i])]
               temp$ub_BPCPmidp[r] <- 1-BPCPmidp[[1]]$lower[which.max(BPCPmidp[[1]]$R>scenarios$tau[i])]
     
               ### CONFIDENCE INTERVALS FOR SIMPLE PROPORTION
               
               ##  analyse those participants with follow up beyond tau 
               nevent <- sum(simdat[which(simdat$followup > scenarios$tau[i] & simdat$time <= scenarios$tau[i]),"event"])
               denom  <- nrow(simdat[which(simdat$followup > scenarios$tau[i]),])
               
               temp$sProp[r] <- nevent/denom
               
               prop          <- binom.confint(x=nevent, n=denom, 
                                              conf.level=1-scenarios$alpha[i], methods=c("exact", "wilson"))
               
               temp$lb_Clopper[r] <- prop[which(prop$method=="exact"),]$lower
               temp$ub_Clopper[r] <- prop[which(prop$method=="exact"),]$upper
               
               temp$lb_W[r]  <- prop[which(prop$method=="wilson"),]$lower
               temp$ub_W[r]  <- prop[which(prop$method=="wilson"),]$upper
               
          }
          
          CPt1 <- scenarios$target_CP[i]
          
          # save information:
          result_i$scenID   <- scenarios$scenID[i]
          result_i$n        <- scenarios$n[i]
          result_i$alpha    <- scenarios$alpha[i]
          result_i$CPt1     <- scenarios$target_CP[i]
          result_i$cProp    <- mean(temp$cProp)
          result_i$sProp    <- mean(temp$sProp)
          
          result_i$method[1]   <- "Greenwood clog-log"
          result_i$lberr[1]    <- mean(temp$lb_GW>CPt1, na.rm=TRUE)
          result_i$uberr[1]    <- mean(temp$ub_GW<CPt1, na.rm=TRUE)
          
          result_i$method[2]   <- "Thomas-Grunkemeier"
          result_i$lberr[2]    <- mean(temp$lb_TG>CPt1, na.rm=TRUE)
          result_i$uberr[2]    <- mean(temp$ub_TG<CPt1, na.rm=TRUE)
          
          result_i$method[3]   <- "BPCP"
          result_i$lberr[3]    <- mean(temp$lb_BPCP>CPt1, na.rm=TRUE)
          result_i$uberr[3]    <- mean(temp$ub_BPCP<CPt1, na.rm=TRUE)
          
          result_i$method[4]   <- "BPCP Mid P"
          result_i$lberr[4]    <- mean(temp$lb_BPCPmidp>CPt1, na.rm=TRUE)
          result_i$uberr[4]    <- mean(temp$ub_BPCPmidp<CPt1, na.rm=TRUE)
          
          result_i$method[5]   <- "Rothman Wilson"
          result_i$lberr[5]    <- mean(temp$lb_RW>CPt1, na.rm=TRUE)
          result_i$uberr[5]    <- mean(temp$ub_RW<CPt1, na.rm=TRUE)
          
          result_i$method[6]   <- "Clopper Pearson"
          result_i$lberr[6]    <- mean(temp$lb_Clopper>CPt1, na.rm=TRUE)
          result_i$uberr[6]    <- mean(temp$ub_Clopper<CPt1, na.rm=TRUE)
          
          result_i$method[7]   <- "Wilson"
          result_i$lberr[7]    <- mean(temp$lb_W>CPt1, na.rm=TRUE)
          result_i$uberr[7]    <- mean(temp$ub_W<CPt1, na.rm=TRUE)
          
          allresults <- rbind(allresults, result_i)
          
     }
     
     return(allresults)

}
```

### Use function to obtain results:

Note: we recommend to use nrep=100,000 but use 1000 here to reduce
computation time.

``` r
example1output <- simulateCIs(n=50, tau=8, cens_max=10, target_CP=0.2, alpha=0.05, nrep=1000)
```

    ## [1] "Scenario 1 of 1; n=50; targetCP=0.2"
    ## [1] "2023-05-30 11:02:51 EDT"

``` r
example1output
```

    ##   scenID  n alpha CPt1    cProp  sProp             method lberr uberr
    ## 1      1 50  0.05  0.2 0.199441 0.2019 Greenwood clog-log 0.029 0.021
    ## 2      1 50  0.05  0.2 0.199441 0.2019 Thomas-Grunkemeier 0.022 0.040
    ## 3      1 50  0.05  0.2 0.199441 0.2019               BPCP 0.009 0.000
    ## 4      1 50  0.05  0.2 0.199441 0.2019         BPCP Mid P 0.016 0.000
    ## 5      1 50  0.05  0.2 0.199441 0.2019     Rothman Wilson 0.023 0.031
    ## 6      1 50  0.05  0.2 0.199441 0.2019    Clopper Pearson 0.007 0.000
    ## 7      1 50  0.05  0.2 0.199441 0.2019             Wilson 0.032 0.000

You can also input multiple values for each parameter:

``` r
example2output <- simulateCIs(n=c(25,50), tau=8, cens_max=10, target_CP=c(0.2, 0.5), alpha=0.05, nrep=1000)
```

    ## [1] "Scenario 1 of 4; n=25; targetCP=0.2"
    ## [1] "2023-05-30 11:03:36 EDT"
    ## [1] "Scenario 2 of 4; n=50; targetCP=0.2"
    ## [1] "2023-05-30 11:04:04 EDT"
    ## [1] "Scenario 3 of 4; n=25; targetCP=0.5"
    ## [1] "2023-05-30 11:04:50 EDT"
    ## [1] "Scenario 4 of 4; n=50; targetCP=0.5"
    ## [1] "2023-05-30 11:05:20 EDT"

``` r
example2output
```

    ##    scenID  n alpha CPt1     cProp  sProp             method lberr uberr
    ## 1       1 25  0.05  0.2 0.2023872 0.2032 Greenwood clog-log 0.028 0.032
    ## 2       1 25  0.05  0.2 0.2023872 0.2032 Thomas-Grunkemeier 0.021 0.071
    ## 3       1 25  0.05  0.2 0.2023872 0.2032               BPCP 0.009 0.000
    ## 4       1 25  0.05  0.2 0.2023872 0.2032         BPCP Mid P 0.018 0.000
    ## 5       1 25  0.05  0.2 0.2023872 0.2032     Rothman Wilson 0.024 0.044
    ## 6       1 25  0.05  0.2 0.2023872 0.2032    Clopper Pearson 0.004 0.000
    ## 7       1 25  0.05  0.2 0.2023872 0.2032             Wilson 0.062 0.000
    ## 8       2 50  0.05  0.2 0.2010829 0.1991 Greenwood clog-log 0.028 0.023
    ## 9       2 50  0.05  0.2 0.2010829 0.1991 Thomas-Grunkemeier 0.022 0.060
    ## 10      2 50  0.05  0.2 0.2010829 0.1991               BPCP 0.010 0.000
    ## 11      2 50  0.05  0.2 0.2010829 0.1991         BPCP Mid P 0.020 0.000
    ## 12      2 50  0.05  0.2 0.2010829 0.1991     Rothman Wilson 0.025 0.034
    ## 13      2 50  0.05  0.2 0.2010829 0.1991    Clopper Pearson 0.009 0.000
    ## 14      2 50  0.05  0.2 0.2010829 0.1991             Wilson 0.030 0.000
    ## 15      3 25  0.05  0.5 0.4934831 0.4974 Greenwood clog-log 0.047 0.025
    ## 16      3 25  0.05  0.5 0.4934831 0.4974 Thomas-Grunkemeier 0.019 0.046
    ## 17      3 25  0.05  0.5 0.4934831 0.4974               BPCP 0.008 0.000
    ## 18      3 25  0.05  0.5 0.4934831 0.4974         BPCP Mid P 0.018 0.001
    ## 19      3 25  0.05  0.5 0.4934831 0.4974     Rothman Wilson 0.037 0.046
    ## 20      3 25  0.05  0.5 0.4934831 0.4974    Clopper Pearson 0.000 0.000
    ## 21      3 25  0.05  0.5 0.4934831 0.4974             Wilson 0.028 0.033
    ## 22      4 50  0.05  0.5 0.5023320 0.4991 Greenwood clog-log 0.034 0.028
    ## 23      4 50  0.05  0.5 0.5023320 0.4991 Thomas-Grunkemeier 0.034 0.037
    ## 24      4 50  0.05  0.5 0.5023320 0.4991               BPCP 0.019 0.005
    ## 25      4 50  0.05  0.5 0.5023320 0.4991         BPCP Mid P 0.029 0.011
    ## 26      4 50  0.05  0.5 0.5023320 0.4991     Rothman Wilson 0.023 0.044
    ## 27      4 50  0.05  0.5 0.5023320 0.4991    Clopper Pearson 0.013 0.016
    ## 28      4 50  0.05  0.5 0.5023320 0.4991             Wilson 0.013 0.016
