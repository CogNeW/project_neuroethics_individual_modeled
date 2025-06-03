#Kant_Processing is a script for preparing data and providing 
#summary measures and plots for intial analysis
rm(list=ls())

source("./r_docs/LoadLibraries.R")
LoadLibraries()

#Import and Transform
files <- list.files(path = "./data/mturk_data_first_50_subjs", pattern = "*.csv", full.names = T)
pilot <- sapply(files, readr::read_csv, simplify=FALSE) 
Mpilot = bind_rows(pilot)

N = 75#Set the min number of trials
n = 54#Set the number of subjects 
cleanPilot = select_(Mpilot,"user_id","condition","total_trials","risk","gain","no_effect", "experimental_treatment_selected", "key_press", "rt")
cleanPilot = na.omit(cleanPilot)
cleanPilot = filter(cleanPilot, total_trials >= N)
cleanPilot = mutate(cleanPilot, ratio = gain / risk)#Gain over Risk
cleanPilot$ratio = round(cleanPilot$ratio,4)
cleanPilot$condition = factor(cleanPilot$condition)


#Check the RT data for outliers
trial_index = data_frame(rep.int(c(1:75), n)) #MAKE SURE THIS IS EQUAL TO THE NUMBER OF SUBJECTS
colnames(trial_index)= c("trial_index")

cleanPilot = cbind(cleanPilot,trial_index)
cleanPilot$trial_index = factor(cleanPilot$trial_index)

cleanPilot = filter(cleanPilot, rt <= 60000 & rt >= 250)

cleanPilot =subset(cleanPilot, with(cleanPilot, user_id %in% names(which(table(user_id)>=N)))) #Remove Subjects who do not have 75 trials

#table(cleanPilot$user_id) Sanity Check
cleanPilot$user_id = factor(cleanPilot$user_id)
RT = select_(cleanPilot,"user_id", "rt", "trial_index")
summary(RT)
##Plots for RT 
library(reshape)
data_long = melt(RT, 
                 id=c("user_id","trial_index"),
                 measure=c("rt"))

summary(data_long)

ggplot(data_long, aes(trial_index, value)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data=mean_se, geom = "pointrange") +
  scale_y_continuous("RT (milliseconds)") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  labs(title="", 
       subtitle="",
       caption="",
       x="Trial Number",
       y="RT (milliseconds)")

datalist = list()
for (i in RT$user_id){
  x = filter(RT, user_id == i)
  dat = summary(x$rt) 
  datalist[[i]] = dat # add it to your list
}

user_id = factor(names(datalist)) #creates vector of user_id's 

for (i in user_id){
  x = filter(RT, user_id == i)
  sum = summary(x)
  print(sum)
  }

rm(data_long)


###Add Survey Data 
survey = read.csv("./data/Qualtrics Data/Risk_Data.csv", header=FALSE, comment.char="#", stringsAsFactors=TRUE)
survey = survey[,-c(1:13)]
questions_strings = survey[1,] #Question ID's 
survey = survey[-1,]
colnames(survey)[1]= 'user_id' 
colnames(survey)[19]= 'Attention_Check' 
survey = filter(survey, Attention_Check == 0)


Demographics = survey[, 97:114]
Demographics = survey[, 97:101]


colnames(Demographics) = c('Sex','Age', 'Race/Ethnicity', 'Education', 'Religion')
Demographics = cbind(Demographics, survey[1])


Cog_Domains = factor(c("Attention (Concentration)", "Motor Function", "Language (spoken)","Long Term Memory","Mood","Self-control","Short Term (Working) Memory"))


cleanPilot = merge(cleanPilot, survey, "user_id")

cleanPilot = merge(cleanPilot, Demographics, "user_id")


Neuroethics_Judgement = cleanPilot
#Check order the levels of the outcome factor to ensure that No to Treatment is  the baseline categorie
Neuroethics_Judgement$experimental_treatment_selected = factor(Neuroethics_Judgement$experimental_treatment_selected, levels = c(FALSE,TRUE), labels = c("No to treatment","Yes to treatment"))

write.csv(Neuroethics_Judgement,'Cleaned_Pilot.csv')

rm(Mpilot, pilot, files, cleanPilot) #Clean up the environment

#How many observatoins does each subject have in each class (yes and no)?
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




TrialMod = by(Neuroethics_Judgement, Neuroethics_Judgement$trial_index, function(x) glm(experimental_treatment_selected ~ risk + gain + ratio, data = x, family = binomial(link = "logit"),control=glm.control(maxit=50)))

sub = filter(Neuroethics_Judgement, user_id == "01rSsifR")


####Logsistic Reg

BaseMod = by(Neuroethics_Judgement, Neuroethics_Judgement$user_id, function(x) glm(experimental_treatment_selected ~ risk + gain, data = x, family = binomial(link = "logit"),control=glm.control(maxit=50)))
Ratio_mod = by(Neuroethics_Judgement, Neuroethics_Judgement$user_id, function(x) glm(experimental_treatment_selected ~ risk + gain + ratio, data = x, family = binomial(link = "logit"),control=glm.control(maxit=50)))
InteractionMod = by(Neuroethics_Judgement, Neuroethics_Judgement$user_id, function(x) glm(experimental_treatment_selected ~ risk + gain + risk*gain, data = x, family = binomial(link = "logit"),control=glm.control(maxit=50)))


user_id = names(BaseMod) #creates vector of user_id's 

#pchisq(84.895-21.849, 73) #Chi squared to see relative differeneces from the null
#anova(mod3$nB8WmX4K, mod4$nB8WmX4K, test="Chisq")

#loop an anove test of the interaction on base and ratio model
datalist = list()
for (i in 1:39){
dat = anova(BaseMod[[i]], Ratio_mod[[i]], test="Chisq")
  datalist[[i]] = dat
}
lapply(datalist, print) 

#drop1(InteractionMod$nB8WmX4K, test="Chisq")
#mod1 = by(Neuroethics_Judgement, Neuroethics_Judgement$user_id, function(x) glm(experimental_treatment_selected ~ gain, data = x, family = binomial(link = "logit"),control=glm.control(maxit=50)))
#mod2 = by(Neuroethics_Judgement, Neuroethics_Judgement$user_id, function(x) glm(experimental_treatment_selected ~ no_effect, data = x, family = binomial(link = "logit"),control=glm.control(maxit=50)))
#mod4 = by(Neuroethics_Judgement, Neuroethics_Judgement$user_id, function(x) glm(experimental_treatment_selected ~ risk*gain, data = x, family = binomial(link = "logit"),control=glm.control(maxit=50)))

#mod3 = by(Neuroethics_Judgement, Neuroethics_Judgement$user_id, function(x) glm(experimental_treatment_selected ~ risk + gain, data = x, family = binomial(link = "logit"),control=glm.control(maxit=50)))



#individual curve fitting for a single predictor
sub = filter(Neuroethics_Judgement, user_id == "nB8WmX4K")
logr <- glm(experimental_treatment_selected ~ risk, data = sub, family = binomial(link = "logit"),control=glm.control(maxit=50))
summary(logr)
ggplot(sub, aes(x=risk, y=experimental_treatment_selected)) + geom_point() + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)


  par(mar = c(4, 4, 1, 1)) # Reduce some of the margins so that the plot fits better
  plot(sub$risk, sub$experimental_treatment_selected)
  curve(predict(logr, data.frame(risk=x), type="response"), add=TRUE) 
  
##########################################
  
#predicted = predict(i, type="response")
#optCutOff = optimalCutoff(ii$experimental_treatment_selected, predicted)[1] 
lapply(BaseMod, summary)
#lapply(mod4[1:10], summary)
#summary(mod3$nB8WmX4K)
#summary(mod4$nB8WmX4K)

#Omnibus =  list(BaseMod, Ratio_mod)
# display results


#lapply(BaseMod,confint) # 95% CI for the coefficients
# exp(coef(m))exponentiated coefficients # 95% CI for exponentiated coefficients exp(confint(m))
#lapply(BaseMod, predict(type="response")) # predicted values
#lapply(BaseMod, residuals( type="deviance")) # residuals



#Extraction of the Coefficients for risk and gain
Coefficients = lapply(BaseMod, coef)
Coefficients = do.call(rbind, Coefficients)
Coefficients = Coefficients[,-1] #Removes intercept
Coefficients = data.frame(Coefficients)
colnames(Coefficients)= c('Risk_Beta', 'Gain_Beta')

Coefficients$user_id = names(BaseMod)

Neuroethics_Judgement= merge(Neuroethics_Judgement, Coefficients, "user_id")
rm(Coefficients)



Coefficients = exp(Coefficients)#exponetiated 
#Rename
XX = exp(Coefficients)#Odds Ratios
colnames(XX) = c('OddR_Risk_Beta', 'OddR_Gain_Beta')
#Round 
XX = round(XX,4)
Coefficients[ , c('Risk_Beta', 'Gain_Beta')] = round(Coefficients[ , c('Risk_Beta', 'Gain_Beta')],4)
#Add Merge parameter
XX$user_id = names(BaseMod)
Coefficients$user_id = names(BaseMod)
Coefficients= merge(Coefficients, XX, "user_id")
#Add to main dataset
Neuroethics_Judgement= merge(Neuroethics_Judgement, Coefficients, "user_id")
rm(Coefficients, XX)


Coefficients = lapply(TrialMod, coef)
Coefficients = do.call(rbind, Coefficients)
Coefficients = Coefficients[,-1] #Removes intercept
Coefficients = data.frame(Coefficients)
colnames(Coefficients)= c('Risk_Beta', 'Gain_Beta', 'Ratio_Beta')

Coefficients$Trial_index = names(TrialMod)

Neuroethics_Judgement= merge(Neuroethics_Judgement, Coefficients, "user_id")
rm(Coefficients)



#Coefficients$Odds_RISK = round(exp(Coefficients[,1]),3)
#Coefficients$Odds_GAIN = round(exp(Coefficients[,2]),3)
#Coefficients = slice(Coefficients, rep(1:n(), each = 75)) #Method only works if all observations are the same length
#we need the exp betas to for interptation purposes
#Add Coefficients to the main CSV
#Neuroethics_Judgement = bind_cols(Neuroethics_Judgement, Coefficients)
names(Neuroethics_Judgement)
#Condition density plots
rm(TEST)
TEST = select_(Neuroethics_Judgement,"user_id","condition","Risk_Beta","Ratio_Beta" ,"Gain_Beta", "Sex", "Age", "Education")
TEST = distinct(Neuroethics_Judgement,user_id, condition, Risk_Beta, Gain_Beta, Ratio_Beta, Ratio_Beta  )


TEST = merge(TEST , Demographics, "user_id")

write.csv(TEST,'39Sub_Betas.csv')

XX = split(TEST, TEST$condition)

#LONG FORMAT OPTION
library(reshape2)
data_long = melt(TEST, id.vars=c("user_id", "condition", "age", "sex" ))

data_long[,"V111"] = as.numeric(data_long[,"V111"])


write.csv(Neuroethics_Judgement,'Cleaned_Pilot.csv')

data_long = melt(Coefficients, id.vars=c("Trial_index"))

View(data_long)

library(lme4)
library(lmerTest)

modTEST <- lmer(value ~ 1 + V110 + V111 + V113  + condition + (1 | condition), data = data_long, REML = FALSE)
summary(modTEST)
##PLOTS##
source("./r_docs/Multiple_plot_function.R")



#one way between anova:
#IV: variable
#IV: condition
#DV: value
t.test(Neuroethics_Judgement, )
aov2 <- aov(lambda ~ sex, data = Neuroethics_Judgement)
summary(aov2)
model.tables(aov2, "means")


#Draft of Main Box Plot Across conditions
ggplot(data_long, aes(variable, y=value, fill=condition)) +
  geom_boxplot()+
  scale_y_continuous("Beta", breaks=seq(-.3, 2,.05), limits=c(-.3, 2)) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  labs(title="", 
       subtitle="an investigation of domain differences",
       caption="Source: CogNew",
       x="Cognitive Domian",
       y="Log Odds")

#Separate Box Plots across conditions
ggplot(data_long, aes(condition, y=value, fill=variable)) +
  geom_boxplot()+
  facet_wrap(~condition, scale="free")


#ID individual
user_id = filter(Neuroethics_Judgement, user_id == "nB8WmX4K")
sub= filter(Neuroethics_Judgement, user_id == "nB8WmX4K")
i = "nB8WmX4K"

##########Omnibus plot 

for (i in user_id){
  sub = filter(Neuroethics_Judgement, user_id == i)
  
  p1 = ggdensity(sub, 
          main = paste("subject: ",toString(i)),
          x = "risk",
          legend.title = "",
          xticks.by = 10,
          add = "mean", rug = TRUE,
          color = "experimental_treatment_selected", fill ="experimental_treatment_selected", palette = c("#f8766d", "#00bfc4"))
  p3 = ggdensity(sub, 
          x = "gain",
          legend = "none",
          xticks.by = 10,
          add = "mean", rug = TRUE,
          color = "experimental_treatment_selected", fill ="experimental_treatment_selected", palette = c("#f8766d", "#00bfc4"))
  p5 = ggdensity(sub, 
          x = "no_effect",
          legend = "none",
          xticks.by = 10,
          add = "mean", rug = TRUE,
          color = "experimental_treatment_selected", fill ="experimental_treatment_selected", palette = c("#f8766d", "#00bfc4"))

  p2=gghistogram(sub, 
            x = "risk",
            ylab = "Number of Occurrences",
            add = "mean", rug = TRUE,
            fill = "experimental_treatment_selected",  palette = c("#f8766d", "#00bfc4"),
            add_density = TRUE, bins = 30,
            legend = "none",
            xticks.by = 10)
  p4=gghistogram(sub, 
            x = "gain",
            ylab = "Number of Occurrences",
            add = "mean", rug = TRUE,
            fill = "experimental_treatment_selected",  palette = c("#f8766d", "#00bfc4"),
            add_density = TRUE, bins = 30,
            legend = "none",
            xticks.by = 10)
  p6=gghistogram(sub, 
            x = "no_effect",
            alpha = .5,
            ylab = "Number of Occurrences",
            add = "mean", rug = TRUE,
            fill = "experimental_treatment_selected",  palette = c("#f8766d", "#00bfc4"),
            add_density = TRUE, bins = 30, 
            legend = "none",
            xticks.by = 10)

  png(filename = paste0(i , ".png"),
    width = 1234, height = 468)
  multiplot(p1, p2, p3, p4, p5, p6, cols=3)
  dev.off()
}

##RATIO PLOTS
for (i in user_id){
  sub = filter(Neuroethics_Judgement, user_id == i)
  
  p1 = ggdensity(sub, 
               x = "ratio",
               legend.title = paste("75 Trials Pilot Subject: ",toString(i)),
               add = "mean", rug = TRUE,
               color = "experimental_treatment_selected", fill ="experimental_treatment_selected", palette = c("#f8766d", "#00bfc4"),
               xticks.by = 1)
  p2 = gghistogram(sub, 
            x = "ratio",
            alpha = .5,
            ylab = "Number of Occurrences",
            add = "mean", rug = TRUE,
            fill = "experimental_treatment_selected",  palette = c("#f8766d", "#00bfc4"),
            add_density = TRUE, bins = 200,
            legend = "none", 
            xticks.by = 1)
  
  png(filename = paste0(i , ".png"),
    width = 1234, height = 468)
  multiplot(p1, p2, cols=1)
  dev.off()
}

#Individual Scatter: subject graph 

for (i in user_id){
  sub = filter(Neuroethics_Judgement, user_id == i)
  
g = ggplot(subject,aes(x = risk , y = gain)) + 
  geom_point(aes(color = experimental_treatment_selected), size = 2.5, alpha = .7) +
  ggtitle(paste("75 Trials Pilot Subject: ",toString(i))) +
  scale_x_continuous("Risk Probabilities", breaks=seq(0,100,5), limits=c(0, 100)) +
  scale_y_continuous("Gain Probabilities", breaks=seq(0,100,5), limits=c(0, 100)) + 
  theme(legend.position="none")
  
  png(filename = paste0(i , ".png"),
      width = 500, height = 372)
  print(g)
  dev.off()
}



ggplot(TEST, aes(Risk_Beta)) +
  geom_density(aes(fill=factor(condition)), alpha=0.8) +
  scale_x_continuous("Beta", breaks=seq(-5,10,1), limits=c(-5, 10))+
  theme(axis.text.x = element_text(angle=65, vjust=0.6))


p2 = ggplot(TEST, aes(Gain_Beta)) +
  geom_density(aes(fill=factor(condition)), alpha=0.8) +
  scale_x_continuous("Beta", breaks=seq(-5,10,1), limits=c(-5, 10))+
  theme(axis.text.x = element_text(angle=65, vjust=0.6))


library(arm)
subjectDignostics = dplyr::filter(Neuroethics_Judgement, user_id == '01rSsifR')

binnedplot(fitted(subjectDignostics), 
           residuals(subjectDignostics, type = "response"), 
           nclass = NULL, 
           xlab = "Expected Values", 
           ylab = "Average residual", 
           main = "Binned residual plot", 
           cex.pts = 0.8, 
           col.pts = 1, 
           col.int = "gray")
