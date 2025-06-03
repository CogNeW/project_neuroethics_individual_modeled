###This just creates a data frame with the pvalues and significance counts for each subject
get_pvalues <- function(data, formula){
  #find the betas
  m.logit <- by(data, data$user_id, 
                function(x) glm(formula = formula, data = x, family = binomial(link = "logit"),control=glm.control(maxit=50)))

#go and get them
m.logit %>% 
  map(broom::tidy) %>% 
  map(pluck("p.value")) -> x 

subjectModel.Pvalue <- as.data.frame(do.call(rbind, x))
subjectModel.Pvalue = tibble::rownames_to_column(subjectModel.Pvalue, "user_id")

subjectModel.Pvalue <- subjectModel.Pvalue %>% 
  mutate(Risk_sig = ifelse(subjectModel.Pvalue$V2 > 0.05, 1, 0),
         Gain_sig = ifelse(subjectModel.Pvalue$V3 > 0.05, 1, 0)) 

subjectModel.Pvalue$Risk_sig <-  factor(subjectModel.Pvalue$Risk_sig, levels = c(0,1), labels = c("Not Significant","Significant"))
subjectModel.Pvalue$Gain_sig <-  factor(subjectModel.Pvalue$Gain_sig, levels = c(0,1), labels = c("Not Significant","Significant"))
return(subjectModel.Pvalue)
}