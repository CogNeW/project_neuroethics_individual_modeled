# Estimate mtDDM parameters using genetic algorithms
#
# Requires R.matlab and GA packages and their dependencies
# R.matlab package is only needed if input data is loaded (as in the example here) from a .mat file
#
# 
# outputs:
  # params: best estimate from grid search. Columns:
  #   [drift slope for attribute 1, drift slope for attribute 2, 
       #   drift latency 1, drift latency 2, barrier/boundary size]
# LL:  log likelihood of the best-fitting grid search
# bic: bayesian information criterion value
# output: details from genetic algorithm search
#
# NOTE!!! input rts need to be in cs. output slopes and latencies are in cs - 
  #     multiply them by 10 to bring into milliseconds
# v 1.1 by njs 	2021-05-12
# v 1.2 change directories to pwd, fix bug in variable assignment storing ga output, add options to ga function that improve search results


remove(list = ls())

# here is where you set parameters that will change from system to system
# directory with data and location to save to
dataDir <- getwd()
saveDir <- getwd()
inputFile <- "exampleData.mat"
outputFile <- "mtDDMoutput.RData"


# below is code for mtDDM estimation

# get necessary packages
library(R.matlab)
library(GA)


# load participant data
# you will need to change this for to load the input file of your choice
data <- readMat(paste(dataDir,inputFile,sep="/"))
values <- data$values
choices <- data$choices
rts <- data$rts
meanrt <- mean(rts)

# constrain estimation for speed
# NOTE: These may need to change depending on input values. Make sure that
# estimations are not hitting upper or lower bounds. You can omit these by
# leaving them empty (e.g. lb = [];) but this will take a very long time to
# estimate.
# order is [drift slope 1, drift slope 2, latency 1, latency 2, bounds (symmetrical)]
# drift rate 1 sets rate of accumulation for values in column one of "values" variable
# drift rate 2 sets rate of accumulation for values in column two of "values" variable
# latency 1 sets entry time for values in column one of "values" variable
# latency 2 sets entry time for values in column two of "values" variable
lbound <- c(0,0,1,1,.5)
ubound <- c(.09,.09,meanrt,meanrt,2)


likefun <- function(params, values, choices, rts) {
  # output = negative likelihood likelihood value for a set of parameters
  
  # bin real RTs (log-spaced, as RT distributions are positively skewed)
  w=hist(log(rts),seq(min(log(rts)),max(log(rts)),length.out=5+1), plot=FALSE)
  RTBins=exp(w$breaks)
  
  possDiffs=unique(values)
  ndiffs = nrow(possDiffs)
  
  # how many simulations of each value difference pair?
    # 1000 is good, 5000 is better, and simulations indicate diminishing 
  # marginal returns after 10000 but 10000 is very slow
  nSims = 1000
  nSteps = round(max(rts)); # max time steps allowed for simulations
  

  rate1 <- params[1]*possDiffs[,1]
  rate2 <- params[2]*possDiffs[,2]
  drates<-cbind(rate1,rate2)
  lat1 <- round(params[3])
  lat2 <- round(params[4])
  ndt <- c(lat1,lat2)
  bound <- params[5]
  rvsNoise = .1
  
  # which latency is first? if equal, this is equivalent to an mDDM
  first <- which.min(ndt)
  second <- which.max(ndt)

  prob <- array(c(rep(0,2), rep(0,length(RTBins)), rep(0, ndiffs)), dim = c(2, length(RTBins), ndiffs))
  for (d in 1:ndiffs) {
    
    # get value signal
    rvs <- rep(0,nSteps)
    rvs[ndt[first]:ndt[second]-1] <- rep(drates[d,first],ndt[second]-ndt[first]+1)
    rvs[ndt[second]:nSteps] <- rep(drates[d,first]+drates[d,second],nSteps-ndt[second]+1)
    
    simChoices <- rep(0,nSims)
    simRTs <- rep(0,nSims)
    for (s in 1:nSims) {
      
      thisrvs <- cumsum(rvs + rnorm(nSteps,0,rvsNoise))

      aboveT <- which(thisrvs>=bound)
      belowB <- which(thisrvs<=-bound)
      
      timeT <- aboveT[1]
      timeB <- belowB[1]
      
      topChoice <- (is.na(timeB) & !is.null(timeT)) |
        (!is.na(timeB) & !is.null(timeT) & timeT < timeB)
      botChoice <- !is.na(timeB) & is.null(timeT) |
        (!is.na(timeB) & !is.null(timeT) & timeT > timeB)
      
      simChoices[s] <- sum(1 * topChoice, -1 * botChoice, na.rm=TRUE)
      simRTs[s] <- sum(timeT * topChoice, timeB * botChoice, na.rm=TRUE)
      
    }
    # which bin do the simulated RTs and choices fall into?
    whichRTBin <- findInterval(simRTs, RTBins)
    if (length(which(whichRTBin > 0) > 0)>0) {
      for (t in unique(whichRTBin[which(whichRTBin > 0)])) {# for experienced RT bins
              prob[1, t, d] <- sum(simChoices==1 & whichRTBin==t)  # top barrier choice
              prob[2, t, d] <- sum(simChoices==-1 & whichRTBin==t)  # bottom barrier choice
      }
    }
  }
  # number -> proportion
  prob <- prob / nSims
  
  
  # evaluate fit to real data
  # find which bin real responses in:
  alltheRespBins <- rep(0,length(choices))
  alltheRespBins[which(choices == 1)] = 1
  alltheRespBins[which(choices == -1)] = 2
  # find which bin real RTs are in:
  allTheRTBins <- findInterval(rts, RTBins)
  # for each real trial, find its likelihood given these params
  like <- rep(0,length(choices))
  for (trial in 1:length(choices)) {
    # find real trial's value bin:
    valBin <- which(possDiffs[,1] == values[trial,1] & possDiffs[,2] == values[trial,2])
    # assign the probability:
    like[trial] <- prob[alltheRespBins[trial], allTheRTBins[trial], valBin]
  }
  LL <- sum(log(like+.Machine$double.eps))
    
  return(LL)
}


# run genetic algorithm optimization (maximizes log likelihood of a set of parameters)
# get some random sets of parameters to include as parents
nSuggest <- 10
suggestedSol <- array(c(rep(0,nSuggest), rep(0, 5)), dim = c(nSuggest, 5))
for (i in 1:5) {
  x<-rnorm(1000,mean(c(ubound[i],lbound[i] )),ubound[i]-lbound[i])
  x <- x[x>lbound[i] & x<ubound[i]]
  suggestedSol[,i] <- sample(x, nSuggest, replace = FALSE, prob = NULL)
}

GA <- ga(type = "real-valued", 
         fitness = function(x) likefun(x,values, choices, rts), 
         lower = c(th = lbound), upper = c(th = ubound))
params <- GA@solution
params[3:4] = round(params[3:4])
LL <- GA@fitnessValue
bic <- -2 * LL + 5 * log(length(choices))

# save the data
save(params, LL, bic, GA, file = paste(saveDir,outputFile,sep="/"))

