#this is another functiont to calcutate the models MAKE SURE YOU HAVE THE CORRECT ONE
risk_lm_test = function(expl) {
  form = as.formula(paste("risk ~", expl))
  lm(form, data = data)
}
gain_lm_test = function(expl) {
  form = as.formula(paste("gain ~", expl))
  lm(form, data = data)
}

plot_fun = function(x, y){
  
  ggplot(data, aes(x=.data[[x]], y=.data[[y]])) +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", se = FALSE) +
    geom_jitter() + 
    labs(x = paste(x), 
         y = paste("Sensitivity: ", stringr::str_sub(y, 1,4)), 
         title  = paste("Sensitivity toward ", stringr::str_sub(y, 1,4)
                        ,"by", x))+
    theme_minimal()+
    stat_regline_equation(label.y = -0.1,
                          aes(label = ..eq.label..))
}