### PERMUTATION HYPOTHESIS TEST ###

# "R datasets package"
library(tidyverse)
library(reshape2)

data %>% 
  dplyr::select(user_id, Risk_Beta, Gain_Beta) %>% 
  distinct() %>% 
  melt(id = c("user_id"), variable.name = "type", value.name = "beta") %>% 
  dplyr::select(type, beta) -> d
# let's add the data into the "data view"
View(d)

# check the names, etc
names(d)
levels(d$type)
# how many observations in each diet?
table(d$type)

# let's look at a boxplot of beta gain by those 2 diets
boxplot(d$beta~d$type, las=1, ylab="beta", 
        xlab="type",main="R & B")

# calculate the difference in sample MEANS
mean(d$beta[d$type=="Risk_Beta"])  # mean for Risk_Beta
mean(d$beta[d$type=="Gain_Beta"])  # mean for Gain_Beta
# lets calculate the absolute diff in means
test.stat1 <- abs(mean(d$beta[d$type=="Risk_Beta"]) - 
                    mean(d$beta[d$type=="Gain_Beta"])) 
test.stat1


# calculate the difference in sample MEDIANS
median(d$beta[d$type=="Risk_Beta"])  # median for Risk_Beta
median(d$beta[d$type=="Gain_Beta"])  # median for Gain_Beta
# lets calculate the absolute diff in medians
test.stat2 <- abs(median(d$beta[d$type=="Risk_Beta"]) - median(d$beta[d$type=="Gain_Beta"]))  #diff in medians
test.stat2

##########################
###  PERMUTATION TEST  ###
##########################

# let's permute...

# for reproducability of results
set.seed(1979)  
# the number of observations to sample
n <- length(d$type)  
# the number of permutation samples to take
P <- 100000
# the variable we will resample from 
#     (note, could use the labels(type) too, and "shuffle this")
variable <- d$beta  

# initialize a matrix to store the permutation data
PermSamples <- matrix(0, nrow=n, ncol=P)
# each column is a permutation sample of data

# now, get those permutation samples, using a loop
# let's take a moment to discuss what that code is doing...
for(i in 1:P){
  PermSamples[,i] <- sample(variable, size= n, replace=FALSE)
}

# we can take a quick look at the first 5 columns of PermSamples
PermSamples[, 1:5]

# now, let's run thru a loop, calculating the test-stats for 
# each sample
### NOTE: could do all these in ONE step...
###       my suggestion, first write code that works, then 
###       tighten it up and make it more efficient
###       (use an 'apply' statement, and write the function 
###        you want)
###    i keep code transparent for teaching purpose

# initialize vectors to store all of the Test-stats:
Perm.test.stat1 <- Perm.test.stat2 <- rep(0, P)

# loop thru, and calculate the test-stats
for (i in 1:P){
  # calculate the perm-test-stat1 and save it
  Perm.test.stat1[i] <- abs(mean(PermSamples[d$type=="Risk_Beta",i]) - 
                               mean(PermSamples[d$type=="Gain_Beta",i]) )
  # calculate the perm-test-stat2 and save it
  Perm.test.stat2[i] <- abs( median(PermSamples[d$type=="Risk_Beta",i]) - 
                               median(PermSamples[d$type=="Gain_Beta",i]) )
}

# before going too far with this, let's remind ourselves of 
# the TEST STATS
test.stat1; test.stat2
# and, take a look at the first 15 permutation-TEST STATS for 1 and 2
round(Perm.test.stat1[1:15], 5)
round(Perm.test.stat2[1:15], 5)

# and, let's calculate the permutation p-value...
# notice how we can ask R a true/false question...(for the first 15)
(Perm.test.stat1 >= test.stat1)[1:15]
# and if we ask for the mean of all of those, it treats 0=FALSE, 1=TRUE
mean((Perm.test.stat1 >= test.stat1)[1:15])

#...calculate the p-value, for all P=100,000
mean(Perm.test.stat1 >= test.stat1)

# and, let's calculate the p-value for 
# option 2 of the test statistic (abs diff in medians)
mean( Perm.test.stat2 >= test.stat2)

# now, remember the difference between "statistical significance" 
# and "scientific significance"
Perm.test.stat1 %>% 
  as.tibble() %>% 
  ggplot(aes(value)) +geom_density(fill = "#ff4d4d") +
  geom_vline(xintercept = test.stat1, linetype="dashed", size=1.5)

Perm.test.stat2 %>% 
  as.tibble() %>% 
  ggplot(aes(value)) +geom_density(fill = "#ff4d4d") +
  geom_vline(xintercept = test.stat2, linetype="dashed", size=1.5)

########################
### SOME EXTRA STUFF ###
########################

# calculate the mean for each group
with(d, tapply(beta, type, mean))
# calculate the difference in means for each group
abs( diff( with(d, tapply(beta, type, mean)) ) )

# and calculate the median for each group
with(d, tapply(beta, type, median))
# and, calculate the difference in medians 
abs( diff( with(d, tapply(beta, type, median)) ) )

## let's take a look at the 3 "Classic" hyp tests we could 
# consider (each of which comes with their own limitations...)

# let's look at the Independent 2-sample t-test
# tests Ho: means are equal
t.test(d$beta~d$type, paired=T, var.eq=F)  
# let's look at the Wilcoxon aka Mann-Whitney U 
# tests Ho: medians are equal
wilcox.test(d$beta~d$type, paired=F)  
# let's look at the Kolmogorov-Smirnov 2-sample test
# tests Ho: distributions are same
ks.test(d$beta[d$type=="Risk_Beta"], d$beta[d$type=="Gain_Beta"], paired=T)     

### produce a plot that shows the sampling distribution, and p-value
### for the Permutation approach, for test.stat1
# let's take a look at a density plot of all the Permutation 
# test-stats, and add in our Observed test stat
plot(density(Perm.test.stat1), 
     xlab=expression( group("|", bar(Yc) - bar(Ym), "|") ) , 
     main="Permutation Test Stats", las=1)
abline(v=test.stat1, col="red", lty="dotted")
text(60,0.0005, "p-value", col="blue", cex=0.7)


plot(density(Perm.test.stat2), 
     xlab=expression( group("|", bar(Yc) - bar(Ym), "|") ) , 
     main="Permutation Test Stats", las=1)
abline(v=test.stat2, col="red", lty="dotted")
text(60,0.0005, "p-value", col="blue", cex=0.7)
