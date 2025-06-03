risk_lm_test = function(expl) {
  form = as.formula(paste("risk_average_marginal_effect ~", expl))
  lm(form, data = data)
}
gain_lm_test = function(expl) {
  form = as.formula(paste("gain_average_marginal_effect ~", expl))
  lm(form, data = data)
}

plot_fun = function(x, y){
  #x = deparse( substitute(x) )
  #y = deparse( substitute(y) )
  
  ggplot(data, aes(x=.data[[x]], y=.data[[y]])) +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", se = FALSE) +
    geom_jitter() + 
    labs(x = paste(x), 
         y = paste("Average Marginal Effect:", stringr::str_sub(y, 1,4)), 
         title  = paste("Average Marginal Effect of", stringr::str_sub(y, 1,4)
 ,"by", x))+
    theme_minimal()+
    stat_regline_equation(label.y = -0.02,
                          aes(label = ..eq.label..))
}
