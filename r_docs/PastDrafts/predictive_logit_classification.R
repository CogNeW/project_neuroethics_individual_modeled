i = 'zrwNx8or'

datalist = list()
for (i in Neuroethics_Judgement$user_id){
x = filter(Neuroethics_Judgement, user_id == i)
dat = table(x$experimental_treatment_selected) #Check Class bias
datalist[[i]] = dat # add it to your list
}
Observation_Counts = do.call(rbind, datalist)
colnames(Observation_Counts) = c('No', 'Yes')
summary(Observation_Counts)
View(Observation_Counts)
rm(i, x,datalist, dat)


iyes = filter(ii, experimental_treatment_selected==TRUE)
ino = filter(ii, experimental_treatment_selected==FALSE)
shapiro.test(ino$risk)
      
libary(ggpubr)
par(mfrow=c(2,2))

ggdensity(ino$risk, 
main = "Density plot of No occurances",
xlab = "risk")

ggdensity(iyes$risk, 
main = "Density plot of Yes occurances",
xlab = "risk")

ggdensity(ino$gain, 
          main = "Density plot of No occurances",
          xlab = "gain")
ggdensity(iyes$gain, 
          main = "Density plot of Yes occurances",
          xlab = "gain")

ggdensity(ino$no_effect, 
          main = "Density plot of No occurances",
          xlab = "no effect")

ggdensity(iyes$no_effect, 
          main = "Density plot of Yes occurances",
          xlab = "no effect")

d <- density(iyes$risk)
plot(d, main="Density")
polygon(d, col="red", border="blue")


ggqqplot(iyes$risk)

ggqqplot(iyes$risk)

#Create Training and Test Samples
      
input_ones = i[which(i$experimental_treatment_selected == 1), ]  # all 1's
input_zeros = i[which(i$experimental_treatment_selected == 0), ]  # all 0's
      
set.seed(1001) #To replicate the sampling
      
      
#0's as 1's equally sampled for training 
      
input_ones_training_rows = sample(1:nrow(input_ones), 0.75*nrow(input_ones)) #  # 1's for training
input_zeros_training_rows = sample(1:nrow(input_zeros), 0.75*nrow(input_ones))  # 0's for training. 
      
# Create Training Data
training_ones = input_ones[input_ones_training_rows, ]  
training_zeros = input_zeros[input_zeros_training_rows, ]
trainingData = rbind(training_ones, training_zeros)  # row bind the 1's and 0's 
      
# Create Test Data
test_ones = input_ones[-input_ones_training_rows, ]
test_zeros = input_zeros[-input_zeros_training_rows, ]
testData = rbind(test_ones, test_zeros)  # row bind the 1's and 0's 
      
LogModel = glm(experimental_treatment_selected ~ risk + gain + ratio, data = trainingData, family=binomial())
m = glm(experimental_treatment_selected ~ risk + gain, data = i, family=binomial())
      
      
summary(LogModel) # display results

predicted = predict(LogModel, testData, type="response")  # predicted scores
predicted = predict(m, i, type="response")  # predicted scores


predicted = predict(subjectModel.2, subject, type="response")  # predicted scores

#just get the r^2 value...

      
#what is a good cut off?
library(InformationValue)
optCutOff = optimalCutoff(testData$experimental_treatment_selected, predicted)[1] 
optCutOff = optimalCutoff(i$experimental_treatment_selected, predicted)[1] 
      
optCutOff = optimalCutoff(subject$experimental_treatment_selected, predicted)[1] 

#library(car) #If was checking for multicollinearity
#vif(m)
      
#checking misclassification, smaller the better
misClassError(i$experimental_treatment_selected, predicted, threshold = optCutOff)
      
#ROC
plotROC(testData$experimental_treatment_selected, predicted)
      
#Concordance is the percentage of pairs, whose scores of actual positive’s are greater than the scores of actual negative’s.
Concordance(testData$experimental_treatment_selected, predicted)
      
#Number of True Pos
sensitivity(testData$experimental_treatment_selected, predicted)
      
#Number of True Neg
specificity(testData$experimental_treatment_selected, predicted)
      
confusionMatrix(testData$experimental_treatment_selected, predicted)
confusionMatrix(ii$experimental_treatment_selected, predicted, threshold = optCutOff)
      

#rm(Mpilot, pilot, files ) #Clean up the environment

#Classification plot
names(ii$experimental_treatment_selected) <- c('Treatment','No Treatment')
ggplot(subject, aes(x = risk , y = gain)) + 
  geom_point(aes(color = experimental_treatment_selected), size = 6, alpha = .5) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  theme(legend.position="none")
