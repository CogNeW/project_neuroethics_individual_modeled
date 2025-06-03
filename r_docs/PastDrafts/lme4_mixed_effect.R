#load the data I send you and this is the solution I found for modeling our data using a logistic mixed model
#scale the risk and gain variables
continuous <- c("risk", "gain")

data[,continuous] <- scale(data[,continuous])  

#intercept
m.base <- glmer(experimental_treatment_selected ~ 1 +
                  (risk + gain | user_id), 
                data = data, family = binomial(link="logit"), 
                control=glmerControl(optimizer="bobyqa"),
                nAGQ = 0)
#add risk
m.1 <- glmer(experimental_treatment_selected ~ risk +
               (risk + gain | user_id), 
             data = data, family = binomial(link="logit"), 
             control=glmerControl(optimizer="bobyqa"),
             nAGQ = 0)
#add gain
m.2 <- glmer(experimental_treatment_selected ~ risk + gain + 
               (risk + gain | user_id), 
             data = data, family = binomial(link="logit"), 
             control=glmerControl(optimizer="bobyqa"),
             nAGQ = 0)
#add condition
m.3 <- glmer(experimental_treatment_selected ~ risk + gain + condition + 
               (risk + gain | user_id), 
             data = data, family = binomial(link="logit"), 
             control=glmerControl(optimizer="bobyqa"),
             nAGQ = 0)
#add interaction risk and gain 
m.4 <- glmer(experimental_treatment_selected ~ risk + gain + condition + risk:gain +
               (risk + gain | user_id), 
             data = data, family = binomial(link="logit"), 
             control=glmerControl(optimizer="bobyqa"),
             nAGQ = 0)
#add interaction risk and condition
m.5 <- glmer(experimental_treatment_selected ~ risk + gain + condition + risk:gain + risk:condition +
               (risk + gain | user_id), 
             data = data, family = binomial(link="logit"), 
             control=glmerControl(optimizer="bobyqa"),
             nAGQ = 0)
#add interaction gain and condition
m.6 <- glmer(experimental_treatment_selected ~ risk + gain + condition + risk:gain + risk:condition + gain:condition +
               (risk + gain | user_id), 
             data = data, family = binomial(link="logit"), 
             control=glmerControl(optimizer="bobyqa"),
             nAGQ = 0)

anova(m.base, m.1, m.2, m.3, m.4, m.5,m.6)





global.estimates <- coef(summary(m.2))
#get our global estimates
estimates <- cbind(global.estimates[1,1], global.estimates[2,1], global.estimates[3,1])
#get the subject random effects
subject.estimates <- add_rownames(ranef(m.4)$user_id[,1:3], var = "user_id")
#add the global estimates to our random deviations giving us subject wise beta estimates
subject.estimates[,2] <- subject.estimates[,2] + global.estimates[1]
subject.estimates[,3] <- subject.estimates[,3] + global.estimates[2]
subject.estimates[,4] <- subject.estimates[,4] + global.estimates[3]
#calculate lambda (loss aversion) lambda > 1 is loss aversive
subject.estimates <- mutate(subject.estimates, lambda = -subject.estimates$risk/subject.estimates$gain)
#round to take a look
subject.estimates[,2:5] <- round(subject.estimates[,2:5], 4)

#give some names and porvide col; this indicates the percent of trials a given subejct said yes...
colnames(subject.estimates) <- c("user_id", "beta.intercept", "beta.risk", "beta.gain", "lambda")
N <- 75
data <- merge(data,subject.estimates, "user_id") %>% 
  mutate(percent = Yes/N*100) 

datasub <- distinct(data, user_id, condition, beta.risk, beta.gain, lambda, percent)


#Now if we take a look at these values, the interpretation tracks what we would expect for a given loss aversion score 
#and the percent of trails where a subject said yes