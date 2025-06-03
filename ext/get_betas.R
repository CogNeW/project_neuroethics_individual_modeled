#this is a draft function for going to get the betas and outputting all the information I want in a clean table
#that has been merged with the data
get_betas <- function(data, formula){
  #find the betas
m.logit <- by(data, data$user_id, 
              function(x) glm(formula = formula, data = x, family = binomial(link = "logit"),control=glm.control(maxit=50)))
#go and get them
map(m.logit, coef) -> subjectModels
subjectModels <- as.data.frame(do.call(rbind, subjectModels))
subjectModels <- tibble::rownames_to_column(subjectModels, "user_id")

other = as.data.frame(do.call(rbind, map(m.logit, broom::glance)))

other = tibble::rownames_to_column(other, "user_id")

user_id = subjectModels$user_id
colin <- list()
#check the variance inflation factors 
for (i in user_id){
  vif <- car::vif(m.logit[[i]])
  colin[[i]] <- vif
}
subjectcolin <- as.data.frame(do.call(rbind, colin))
subjectcolin <- tibble::rownames_to_column(subjectcolin, "user_id")
#remove subject if they have vif for risk or gain above 5
#this is the part I need to fix, it is too specific to the data in this example
subjectcolin <- subjectcolin %>% 
  mutate(index = ifelse(subjectcolin[,2] > 5 | subjectcolin[,3] > 5, 1, 0)) %>% 
  filter(index == 0)
#merge the subject models extracted above and the index
subjectModels <- merge(subjectModels, subjectcolin, "user_id")
subjectModels <- merge(subjectModels, other, "user_id")

#create data that only contains subjects after VIF data exculsion
data <- merge(data, subjectModels, "user_id")
return(data)
} 