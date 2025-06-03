# via https://stat.ethz.ch/pipermail/r-sig-mixed-models/2008q2/000874.html
# Calculate the median odds ratio
# -------------------------------
MOR <- function(my.var, digits = 2)
          { # MOR arguments: my.var = variance associated with level 2 clustering variable
            # digits = number of decimal places to which MOR value will be rounded.

            Median.OR <- round(exp(sqrt(2*my.var)*qnorm(.75)), digits)
            paste("Median Odds-Ratio (MOR) = ", Median.OR)  }

IOR <- function(my.var, my.coef, my.diff, digits = 2)
          { # IOR arguments: my.var = variance associated with level 2 clustering variable
            # my.coef = fixed effect coefficient associated with a level 2 covariate
            # my.diff = difference between 2 different possible values on the level
            #           2 covariate associated with my.coef
            # digits = number of decimal places to which MOR value will be rounded.

            IOR.LL <- round(exp(my.coef*my.diff + sqrt(2*my.var)*qnorm(.10)), digits)
            IOR.UL <- round(exp(my.coef*my.diff + sqrt(2*my.var)*qnorm(.90)), digits)
            print("80% Interval Odds Ratio lower & upper limits (IOR.LL & IOR.UL)")
            cbind(my.var, my.coef, my.diff, IOR.LL, IOR.UL)
            }