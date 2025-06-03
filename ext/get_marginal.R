#this is a draft function for going to get the betas and outputting all the information I want in a clean table
#that has been merged with the data
get_marginal <- function(data, formula){
  #find the marginal effect
  m.marginal = by(logit.data, logit.data$user_id, function(x) logitmfx(treatment ~ risk.c + gain.c, control=glm.control(maxit=50), data = x))
  user_id = names(m.marginal)
  marginal.effect <- list()
  
  for (i in user_id){
    risk = list(m.marginal[[i]]$mfxest[1,])
    gain = list(m.marginal[[i]]$mfxest[2,])
    marginal.effect[[i]] = cbind(risk, gain)
  }
  marginal.effect.risk = map_df(marginal.effect, pluck(1))
  marginal.effect.gain = map_df(marginal.effect, pluck(2))
  
  subjectModels.marginal <- cbind.data.frame(user_id, marginal.effect.risk, marginal.effect.gain)  
  colnames(subjectModels.marginal) <- c("user_id", "risk_average_marginal_effect", "risk_Std_Err", "risk_z", "risk_P>|z|", "gain_average_marginal_effect", "gain_Std_Err", 
                                        "gain_z", "gain_P>|z|")

  #create data that only contains subjects after VIF data exculsion
  data <- merge(data, subjectModels.marginal, "user_id")
  return(data)
}
