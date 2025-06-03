# kim_analysis_stddm - conducts parameter estimation and model comparisons for 
# neuroethics project
#
# Copyright (C) 2024 Kianté Fernandez, <kiantefernan@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Record of Revisions
#
# Date            Programmers                         Descriptions of Change
# ====         ================                       ======================
# 2024/08/20      Kianté Fernandez                       wrote code v1

# Load required libraries
library(runjags) # Interface Utilities, Model Templates, Parallel Computing Methods and Additional Distributions for MCMC Models in JAGS, CRAN v2.2.0-3
library(dplyr) # A Grammar of Data Manipulation, CRAN v1.0.9
library(tidyr) # Tidy Messy Data, CRAN v1.2.0
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics, CRAN v3.4.4
library(tidybayes) # Tidy Data and 'Geoms' for Bayesian Models, CRAN v3.0.2
library(ggdist) # Visualizations of Distributions and Uncertainty, CRAN v3.0.1
library(loo) # Efficient Leave-One-Out Cross-Validation and WAIC for Bayesian Models, CRAN v2.5.1
library(GGally) # Extension to 'ggplot2', CRAN v2.1.2
library(RWiener) # Wiener Process Distribution Functions, CRAN v1.3-3
library(tidyverse) # Easily Install and Load the 'Tidyverse', CRAN v1.3.1
library(patchwork) # The Composer of Plots, CRAN v1.1.1

# Load and prepare data
behavior_path <- "C:/Users/Brian Kim/Documents/OneDrive - Drexel University/Documents/classes/cognew/Neuroethics/IndividualDifferenceRiskBenefit/DATA_ANALYSIS/data/FINAL_RT_EXCLUDED_DATA_20241014.RData"
load(behavior_path)

# Data preprocessing
risk_data <- risk_data %>%
  mutate(
    rt = rt / 1000, # Convert reaction time to seconds
    response = ifelse(experimental_treatment_selected == "Yes to treatment", 1, 0),
    RT = if_else(response == 0, -rt, rt), # Negative RT for 'No' responses
    risk = risk * 0.01, # Convert percentages to proportions
    gain = gain * 0.01,
    no_effect = no_effect * 0.01,
    idxP = as.numeric(as.factor(user_id)) # Create numeric subject index
  ) 
# %>%
  # filter(idxP %in% 1:5) # Select first 5 subjects (change this for running lots of people)

# Prepare data for the JAGS model
dat <- with(
  risk_data,
  dump.format(list(N = length(RT), y = RT, idxP = idxP, gain = gain, risk = risk, rt = rt, ns = length(unique(idxP))))
)

# Define initial values for the MCMC chains
inits <- list(
  dump.format(list(
    alpha.mu = 2.4, alpha.pr = 0.05, theta.mu = 0.15, theta.pr = 0.05,
    b_gain.mu = 0.1, b_gain.pr = 0.05, b_risk.mu = 0.05, b_risk.pr = 0.05, bias.mu = 0.3, bias.kappa = 1,
    .RNG.name = "base::Mersenne-Twister", .RNG.seed = 6666
  )),
  dump.format(list(
    alpha.mu = 2.2, alpha.pr = 0.05, theta.mu = 0.01, theta.pr = 0.05,
    b_gain.mu = 0.3, b_gain.pr = 0.05, b_risk.mu = 0.1, b_risk.pr = 0.05, bias.mu = 0.5, bias.kappa = 1,
    .RNG.name = "base::Wichmann-Hill", .RNG.seed = 1234
  )),
  dump.format(list(
    alpha.mu = 2, alpha.pr = 0.05, theta.mu = 0.1, theta.pr = 0.05,
    b_gain.mu = 0.3, b_gain.pr = 0.05, b_risk.mu = 0.01, b_risk.pr = 0.05, bias.mu = 0.4, bias.kappa = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 99999
  ))
)

# Define parameters to monitor (for both models)
monitor_0 <- c("b_gain.mu", "b_risk.mu", "alpha.mu", "theta.mu", "bias.mu", "deviance", "ld_comp")
monitor_1 <- c("b_gain.mu", "b_risk.mu", "alpha.mu", "theta.mu", "bias.mu", "time.mu", "deviance", "ld_comp")

# Function to run JAGS model
run_jags_model <- function(model_file, model_name, monitor) {
  results <- run.jags(
    model = model_file,
    monitor = monitor,
    data = dat,
    n.chains = 3,
    inits = inits,
    plots = FALSE,
    method = "parallel",
    module = "wiener",
    burnin = 1000,
    sample = 5000,
    thin = 1
  )

  mcmc_samples <- coda::as.mcmc.list(results)
  posterior_samples <- mcmc_samples %>% tidybayes::tidy_draws()

  list(
    results = results,
    mcmc_samples = mcmc_samples,
    posterior_samples = posterior_samples
  )
}

# Run models (change the paths for the files I sent)
m0_results <- run_jags_model("C:/Users/Brian Kim/Documents/OneDrive - Drexel University/Documents/classes/cognew/Neuroethics/IndividualDifferenceRiskBenefit/DATA_ANALYSIS/analysis/finalDDM/stddm_model_M0.txt", "M0", monitor_0)
m1_results <- run_jags_model("C:/Users/Brian Kim/Documents/OneDrive - Drexel University/Documents/classes/cognew/Neuroethics/IndividualDifferenceRiskBenefit/DATA_ANALYSIS/analysis/finalDDM/stddm_model_M1.txt", "M1", monitor_1)

# model comparison indices
calculate_indices <- function(posterior_samples, mcmc_samples) {
  # Calculate DIC (sigh... just do it by hand I guess)
  deviance_samples <- posterior_samples$deviance
  D_bar <- mean(deviance_samples)
  D_hat <- min(deviance_samples)
  pD <- D_bar - D_hat
  DIC <- D_hat + 2 * pD

  # Calculate LOO and WAIC
  log_lik_list <- lapply(mcmc_samples, function(chain) {
    matrix(chain[, grep("ld_comp", colnames(chain))], ncol = nrow(risk_data))
  })
  log_lik <- do.call(rbind, log_lik_list)
  chain_id <- rep(seq_along(mcmc_samples), each = nrow(mcmc_samples[[1]]))
  r_eff <- loo::relative_eff(exp(log_lik), chain_id = chain_id, cores = 4)
  loo_result <- loo::loo(log_lik, r_eff = r_eff, cores = 4)
  waic_result <- loo::waic(log_lik)

  list(DIC = DIC, LOO = loo_result, WAIC = waic_result)
}

# Compute indices for models
m0_indices <- calculate_indices(m0_results$posterior_samples, m0_results$mcmc_samples)
m1_indices <- calculate_indices(m1_results$posterior_samples, m1_results$mcmc_samples)

# Print model comparison results
cat("\nModel Comparison Results:\n")
cat("-------------------------\n")
cat("M0 (without time effect):\n")
cat(sprintf("DIC:  %.2f\n", m0_indices$DIC))
cat(sprintf(
  "WAIC: %.2f (SE: %.2f)\n", m0_indices$WAIC$estimates["waic", "Estimate"],
  m0_indices$WAIC$estimates["waic", "SE"]
))
cat(sprintf(
  "LOO:  %.2f (SE: %.2f)\n", m0_indices$LOO$estimates["looic", "Estimate"],
  m0_indices$LOO$estimates["looic", "SE"]
))

cat("\nM1 (with time effect):\n")
cat(sprintf("DIC:  %.2f\n", m1_indices$DIC))
cat(sprintf(
  "WAIC: %.2f (SE: %.2f)\n", m1_indices$WAIC$estimates["waic", "Estimate"],
  m1_indices$WAIC$estimates["waic", "SE"]
))
cat(sprintf(
  "LOO:  %.2f (SE: %.2f)\n", m1_indices$LOO$estimates["looic", "Estimate"],
  m1_indices$LOO$estimates["looic", "SE"]
))

# Compare models using LOO function (model m1 stddm should win here)
loo_comparison <- loo::loo_compare(m0_indices$LOO, m1_indices$LOO)
print(loo_comparison)

browser()

# NOTE: Beyond here is experimental and you can/should play around with plotting
# Create comparison plots
create_comparison_plot <- function(m0_samples, m1_samples, params) {
  m0_data <- m0_samples %>%
    select(all_of(params)) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(model = "M0")

  m1_data <- m1_samples %>%
    select(all_of(params)) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(model = "M1")

  combined_data <- bind_rows(m0_data, m1_data)

  ggplot(combined_data, aes(x = value, y = parameter, fill = model)) +
    stat_halfeye(alpha = 0.5, position = position_dodge(width = 0.5)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = "Comparison of Posterior Distributions", x = "Parameter Value", y = "") +
    theme_tidybayes() +
    scale_fill_manual(values = c("M0" = "skyblue", "M1" = "lightgreen"))
}

params_to_compare <- c("b_gain.mu", "b_risk.mu", "alpha.mu", "theta.mu", "bias.mu")
comparison_plot <- create_comparison_plot(m0_posterior, m1_posterior, params_to_compare)
print(comparison_plot)
ggsave("Figures/ComparisonPlot.png", plot = comparison_plot, dpi = 300, width = 8, height = 6, units = "in")


# Generate predictions and create PPC plots
# TODO: we need a proper function todo the simulations for the stddm I just use DDM simulator for now

library(rtdist)

generate_predictions_standard_ddm <- function(posterior_samples, data, n_pred = 100) {
  n_samples <- nrow(posterior_samples)
  n_obs <- nrow(data)

  sample_indices <- sample(1:n_samples, n_pred, replace = TRUE)
  predictions <- matrix(NA, nrow = n_obs, ncol = n_pred)

  for (i in 1:n_pred) {
    sample <- posterior_samples[sample_indices[i], ]
    for (j in 1:n_obs) {
      alpha <- sample$alpha.mu
      tau <- sample$theta.mu
      beta <- sample$bias.mu
      delta <- sample$b_gain.mu * data$gain[j] + sample$b_risk.mu * data$risk[j]

      pred <- rwiener(1, alpha, tau, beta, delta)
      predictions[j, i] <- ifelse(pred$resp == "upper", pred$q, -pred$q)
    }
  }

  return(predictions)
}

generate_predictions_time_ddm <- function(posterior_samples, data, n_pred = 100) {
 
  # Initialize CPP Files
  # Install required packages if necessary:
  want = c("DEoptim", "Rcpp", "plyr", "parallel", "RcppParallel","stats4","pracma")
  have = want %in% rownames(installed.packages())
  if ( any(!have) ) { install.packages( want[!have] ) }
  # Now load them all
  suppressPackageStartupMessages(lapply(want, require, character.only = TRUE))
  RcppParallel::setThreadOptions(numThreads = 1) #this is critical for mclapply
  
  #load the c++ file containing the functions to simulate the delayed DDM
  # will give some errors that can be ignored
  sourceCpp(file.path(getwd(), "2ddm_r_cpp_M1.cpp"))
  
  # Set up parameters for Stimulation
  n_samples <- nrow(posterior_samples)
  n_obs <- nrow(data)
  
  sample_indices <- sample(1:n_samples, n_pred, replace = TRUE)
  predictions <- matrix(NA, nrow = n_obs, ncol = n_pred)
  
  for (i in 1:n_pred) {
    print(sprintf("i: %d", i))
    sample <- posterior_samples[sample_indices[i], ]
    for (j in 1:n_obs) {
      alpha <- sample$alpha.mu
      tau <- sample$theta.mu
      beta <- sample$bias.mu
      delta <- sample$b_gain.mu * data$gain[j] + sample$b_risk.mu * data$risk[j]
      
      #ddm2_parallel(d_v, d_h, thres, nDT, tIn, bias, vd, hd, 1, N)
      preds = ddm2_parallel(sample$b_gain.mu,sample$b_risk.mu,sample$alpha.mu,sample$theta.mu,sample$time.mu,sample$bias.mu,data$gain[j],data$risk[j],1,1)
      predictions[j, i] <- preds
    }
  }
  
  return(predictions)
}


### -----
#load the c++ file containing the functions to simulate the delayed DDM
# will give some errors that can be ignored
sourceCpp(file.path(main_folder, "2ddm_r_cpp_M1.cpp"))

# Need to load the correct data

DDataSim<- data.frame(hd = risk_data$gain, td = risk_data$risk, subject = risk_data$idxP )
### Simulations ###

# parameters to recover with the fitting

DataSim<-DDataSim
subj <- unique(DataSim$subject)
countSub=1
for (s in subj){
  
  wt = subset(results, X == paste0("b2.p[", toString(s), "]"))$Mean
  wh = subset(results, X == paste0("b1.p[", toString(s), "]"))$Mean
  boundary = subset(results, X == paste0("alpha.p[", toString(s), "]"))$Mean
  nDT = subset(results, X == paste0("theta.p[", toString(s), "]"))$Mean
  tHin = subset(results, X == paste0("time[", toString(s), "]"))$Mean
  bias = subset(results, X == paste0("bias[", toString(s), "]"))$Mean
  
  DataSim$wt[DataSim$subject==s] = wt
  DataSim$wh[DataSim$subject==s] = wh
  DataSim$boundary[DataSim$subject==s] = boundary
  DataSim$nDT[DataSim$subject==s] = nDT
  DataSim$tHin[DataSim$subject==s] = tHin
  DataSim$bias[DataSim$subject==s] = bias
  
  vd <- (DataSim$td[DataSim$subject==s])
  hd <- (DataSim$hd[DataSim$subject==s])
  
  for (i in 1:length(vd)) {
    #ddm2_parallel(d_v, d_h, thres, nDT, tIn, bias, vd, hd, 1, N)
    rts = ddm2_parallel(wt,wh,boundary,nDT,tHin,bias,vd[i],hd[i],1,1)
    if (rts<0){
      DataSim$choice[DataSim$subject==s & DataSim$td ==vd[i] & DataSim$hd ==hd[i]] <- 0
      DataSim$rt[DataSim$subject==s & DataSim$td ==vd[i] & DataSim$hd ==hd[i]] <- -rts
    }else{
      DataSim$choice[DataSim$subject==s & DataSim$td ==vd[i] & DataSim$hd ==hd[i]] <- 1
      DataSim$rt[DataSim$subject==s & DataSim$td ==vd[i] & DataSim$hd ==hd[i]] <- rts
    }
  }
  countSub=countSub+1
}

create_ppc_plot <- function(predictions, observed_data, model_name) {
  predicted_data <- tibble(
    RT = as.vector(predictions),
    draw = rep(1:ncol(predictions), each = nrow(predictions))
  )

  ggplot() +
    stat_density(
      data = predicted_data, aes(x = abs(RT), group = draw),
      geom = "line", position = "identity", alpha = 0.1, color = "lightblue"
    ) +
    stat_density(
      data = predicted_data, aes(x = abs(RT)),
      geom = "line", position = "identity", alpha = 1, color = "blue", linewidth = 1.5
    ) +
    stat_density(
      data = observed_data, aes(x = abs(RT)),
      geom = "line", position = "identity", alpha = 1, color = "red", linewidth = 1.5
    ) +
    labs(
      title = paste("Posterior Predictive Check:", model_name),
      x = "Response Time (seconds)",
      y = "Density"
    ) +
    theme_tidybayes() +
    scale_x_continuous(limits = c(0,7.5)) + 
    # scale_x_continuous(limits = c(0, quantile(abs(c(predicted_data$RT, observed_data$RT)), 0.99))) +
    annotate("text",
      x = Inf, y = Inf, label = "Observed", color = "red",
      hjust = 1.1, vjust = 2, size = 4
    ) +
    annotate("text",
      x = Inf, y = Inf, label = "Predicted (average)", color = "blue",
      hjust = 1.1, vjust = 3.5, size = 4
    ) +
    annotate("text",
      x = Inf, y = Inf, label = "Predicted (samples)", color = "lightblue",
      hjust = 1.1, vjust = 5, size = 4
    )
}

# Generate predictions
set.seed(6969)
m0_predictions <- generate_predictions_standard_ddm(m0_posterior, risk_data, n_pred = 50)
m1_predictions <- generate_predictions_time_ddm(m1_posterior, risk_data, n_pred = 50)

# Create PPC plots
observed_data <- risk_data %>% select(RT)
m0_ppc_plot <- create_ppc_plot(m0_predictions, observed_data, "M0 (without time effect)")
m1_ppc_plot <- create_ppc_plot(m1_predictions, observed_data, "M1 (with time effect)")

print(m0_ppc_plot)
print(m1_ppc_plot)

# combined?
combined_ppc_plot <- ggplot() +
  stat_density(
    data = tibble(RT = as.vector(m0_predictions), model = "M0"),
    aes(x = abs(RT), color = model), geom = "line", alpha = 0.5
  ) +
  stat_density(
    data = tibble(RT = as.vector(m1_predictions), model = "M1"),
    aes(x = abs(RT), color = model), geom = "line", alpha = 0.5
  ) +
  stat_density(
    data = observed_data, aes(x = abs(RT)),
    geom = "line", color = "red", linewidth = 1.5
  ) +
  labs(
    title = "Posterior Predictive Check",
    x = "Response Time (seconds)",
    y = "Density"
  ) +
  theme_tidybayes() +
  scale_color_manual(values = c("M0" = "skyblue", "M1" = "lightgreen", "Observed" = "red")) +
  scale_x_continuous(limits = c(0, quantile(abs(c(m0_predictions, m1_predictions, observed_data$RT)), 0.99))) +
  theme(legend.position = "bottom")

print(combined_ppc_plot)

# Calculate some summary statistics
calculate_summary_stats <- function(predictions, observed_data) {
  observed_summary <- observed_data %>%
    summarise(
      mean_rt = mean(abs(RT)),
      median_rt = median(abs(RT)),
      sd_rt = sd(abs(RT))
    )

  predicted_summary <- tibble(RT = as.vector(predictions)) %>%
    group_by(draw = rep(1:ncol(predictions), each = nrow(predictions))) %>%
    summarise(
      mean_rt = mean(abs(RT)),
      median_rt = median(abs(RT)),
      sd_rt = sd(abs(RT))
    ) %>%
    summarise(across(everything(), list(mean = mean, lower = ~ quantile(., 0.025), upper = ~ quantile(., 0.975))))

  list(observed = observed_summary, predicted = predicted_summary)
}

m0_summary <- calculate_summary_stats(m0_predictions, observed_data)
m1_summary <- calculate_summary_stats(m1_predictions, observed_data)

# Print summary statistics
cat("\nObserved RT Summary:\n")
print(m0_summary$observed)

cat("\nM0 Predicted RT Summary (mean and 95% CI):\n")
print(m0_summary$predicted)
cat("\nM1 Predicted RT Summary (mean and 95% CI):\n")
print(m1_summary$predicted)

# Create a comparison plot with overlaid data
create_comparison_plot_rt <- function(observed_data, predicted_data, x_var, y_var, x_label, y_label, y_limits = NULL, convert_predicted_to_seconds = FALSE) {
  observed_summary <- observed_data %>%
    mutate({{ x_var }} := {{ x_var }} * 100) %>%
    group_by({{ x_var }}, experimental_treatment_selected) %>%
    summarise(
      m_left = mean({{ y_var }}),
      se = sqrt(var({{ y_var }}) / n()),
      .groups = "drop"
    )
  
  predicted_summary <- predicted_data %>%
    mutate({{ x_var }} := {{ x_var }} * 100) %>%
    group_by({{ x_var }}, experimental_treatment_selected) %>%
    summarise(
      m_left = if (convert_predicted_to_seconds) mean({{ y_var }}) / 1000 else mean({{ y_var }}),
      se = if (convert_predicted_to_seconds) sqrt(var({{ y_var }})) / (1000 * sqrt(n())) else sqrt(var({{ y_var }}) / n()),
      .groups = "drop"
    )

  plot <- ggplot() +
    geom_pointrange(
      data = observed_summary,
      aes(x = {{ x_var }}, y = m_left, ymin = m_left - se, ymax = m_left + se, color = "Observed", shape = experimental_treatment_selected)
    ) +
    geom_line(data = observed_summary, aes(x = {{ x_var }}, y = m_left, color = "Observed", shape = experimental_treatment_selected, linetype = experimental_treatment_selected), size = 0.6) +
    geom_pointrange(
      data = predicted_summary,
      aes(x = {{ x_var }}, y = m_left, ymin = m_left - se, ymax = m_left + se, color = "Predicted", shape = experimental_treatment_selected)
    ) +
    geom_line(data = predicted_summary, aes(x = {{ x_var }}, y = m_left, color = "Predicted", shape = experimental_treatment_selected, linetype = experimental_treatment_selected), size = 0.6) +
    theme_classic() +
    geom_smooth(data=observed_summary, aes(x = {{ x_var }}, y = m_left), method="lm") +
    geom_smooth(data=predicted_summary, aes(x = {{ x_var }}, y = m_left), method="lm") +
    # geom_hline(yintercept = .5, linetype = "dashed", size = .2) +
    labs(y = y_label, x = x_label, color = "Data") +
    scale_color_manual(values = c("Observed" = "blue", "Predicted" = "red")) +
    theme(text = element_text(size = 15), legend.position = "bottom")

  if (!is.null(y_limits)) {
    plot <- plot + scale_y_continuous(limits = y_limits)
  }

  return(plot)
}

# Create a comparison plot with overlaid data
create_comparison_plot_choice <- function(observed_data, predicted_data, x_var, y_var, x_label, y_label, y_limits = NULL, convert_predicted_to_seconds = FALSE) {
  observed_summary <- observed_data %>%
    mutate({{ x_var }} := {{ x_var }} * 100) %>%
    group_by({{ x_var }}) %>%
    summarise(
      m_left = mean({{ y_var }}),
      se = sqrt(var({{ y_var }}) / n()),
      .groups = "drop"
    )
  
  predicted_summary <- predicted_data %>%
    mutate({{ x_var }} := {{ x_var }} * 100) %>%
    group_by({{ x_var }}) %>%
    summarise(
      m_left = if (convert_predicted_to_seconds) mean({{ y_var }}) / 1000 else mean({{ y_var }}),
      se = if (convert_predicted_to_seconds) sqrt(var({{ y_var }})) / (1000 * sqrt(n())) else sqrt(var({{ y_var }}) / n()),
      .groups = "drop"
    )

  plot <- ggplot() +
    geom_pointrange(
      data = observed_summary,
      aes(x = {{ x_var }}, y = m_left, ymin = m_left - se, ymax = m_left + se, color = "Observed")
    ) +
    geom_line(data = observed_summary, aes(x = {{ x_var }}, y = m_left, color = "Observed"), size = 0.6) +
    geom_pointrange(
      data = predicted_summary,
      aes(x = {{ x_var }}, y = m_left, ymin = m_left - se, ymax = m_left + se, color = "Predicted")
    ) +
    geom_line(data = predicted_summary, aes(x = {{ x_var }}, y = m_left, color = "Predicted"), size = 0.6) +
    theme_classic() +
    geom_smooth(data=observed_summary, aes(x = {{ x_var }}, y = m_left), method="lm") +
    geom_smooth(data=predicted_summary, aes(x = {{ x_var }}, y = m_left), method="lm") +
    # geom_hline(yintercept = .5, linetype = "dashed", size = .2) +
    labs(y = y_label, x = x_label, color = "Data") +
    scale_color_manual(values = c("Observed" = "blue", "Predicted" = "red")) +
    theme(text = element_text(size = 15), legend.position = "bottom")
  
  if (!is.null(y_limits)) {
    plot <- plot + scale_y_continuous(limits = y_limits)
  }
  
  return(plot)
}

# create fake simulated predictions to mess around with plotting for now.

# predicted_data <- risk_data %>%
#   mutate(
#     response = rbinom(n(), 1, plogis(0.5 + 0.02 * risk - 0.02 * gain)),
#     rt = rnorm(n(), mean = 1500 + 5 * risk - 5 * gain, sd = 200)
#   )

# Calculate the row averages of m0_predictions
row_averages <- rowMeans(m0_predictions)

# Create the new m0_choice_predictions table
m0_choice_predictions <- risk_data %>%
  mutate(
    response = ifelse(row_averages >= 0, 1, 0),  # 1 for positive, 0 for negative
    rt = abs(row_averages)  # Absolute value of the row averages
  )

# Calculate the row averages of m0_predictions
row_averages <- rowMeans(m1_predictions)

# Create the new m0_choice_predictions table
m1_choice_predictions <- risk_data %>%
  mutate(
    response = ifelse(row_averages >= 0, 1, 0),  # 1 for positive, 0 for negative
    rt = abs(row_averages)  # Absolute value of the row averages
  )


# Create plots
m0_choice_risk_plot <- create_comparison_plot_choice(
  risk_data,
  m0_choice_predictions,
  risk,
  response,
  "Probability of Worsening (%)",
  "Probability of Choosing Experimental Treatment",
  y_limits = c(0, 1)
)

m0_choice_gain_plot <- create_comparison_plot_choice(
  risk_data,
  m0_choice_predictions,
  gain,
  response,
  "Probability of Curing (%)",
  "Probability of Choosing Experimental Treatment",
  y_limits = c(0, 1)
)

m0_rt_risk_plot <- create_comparison_plot_rt(
  risk_data,
  m0_choice_predictions,
  risk,
  rt,
  "Probability of Worsening (%)",
  "Response Time (seconds)",
  convert_predicted_to_seconds = FALSE
)

m0_rt_gain_plot <- create_comparison_plot_rt(
  risk_data,
  m0_choice_predictions,
  gain,
  rt,
  "Probability of Curing (%)",
  "Response Time (seconds)",
  convert_predicted_to_seconds = FALSE
)

# Create plots
m1_choice_risk_plot <- create_comparison_plot_choice(
  risk_data,
  m1_choice_predictions,
  risk,
  response,
  "Probability of Worsening (%)",
  "Probability of Choosing Experimental Treatment",
  y_limits = c(0, 1)
)

m1_choice_gain_plot <- create_comparison_plot_choice(
  risk_data,
  m1_choice_predictions,
  gain,
  response,
  "Probability of Curing (%)",
  "Probability of Choosing Experimental Treatment",
  y_limits = c(0, 1)
)

m1_rt_risk_plot <- create_comparison_plot_rt(
  risk_data,
  m1_choice_predictions,
  risk,
  rt,
  "Probability of Worsening (%)",
  "Response Time (seconds)",
  convert_predicted_to_seconds = FALSE
)

m1_rt_gain_plot <- create_comparison_plot_rt(
  risk_data,
  m1_choice_predictions,
  gain,
  rt,
  "Probability of Curing (%)",
  "Response Time (seconds)",
  convert_predicted_to_seconds = FALSE
)

m1_rt_gain_plot_smoothed = m1_rt_gain_plot 
print(m1_rt_gain_plot_smoothed)

print(m0_choice_risk_plot)
print(m0_choice_gain_plot)
print(m0_rt_risk_plot)
print(m0_rt_gain_plot)

print(m1_choice_risk_plot)
print(m1_choice_gain_plot)
print(m1_rt_risk_plot)
print(m1_rt_gain_plot)

# Combine all plots into a panel figure
combined_plot <- (choice_risk_plot + choice_gain_plot) / (rt_risk_plot + rt_gain_plot) +
  plot_annotation(
    title = "Model Predictions vs. Observed Data",
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5))
  ) &
  theme(legend.position = "bottom")

print(combined_plot)

# Extract DDM Metrics

library(dplyr)
library(stringr)

# Step 1: Identify columns with square brackets in m1_posterior
# Extract columns that match the pattern with ".p[XX]" or just "[XX]"
selected_columns <- m1_posterior %>%
  select(matches("\\[\\d+\\]$"))

# Step 2: Calculate the averages for each column
averages <- selected_columns %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# Step 3: Extract both metric_name and participant_id, handling both ".p[XX]" and "[XX]"
column_info <- tibble(
  full_colname = colnames(averages),
  participant_id = str_extract(colnames(averages), "(?<=\\[)\\d+(?=\\])"),  # Extract participant ID if present in square brackets
  metric_name = str_replace(colnames(averages), "\\.?p?\\[\\d+\\]$", "")  # Remove either ".p[XX]" or "[XX]"
)

# Pivot the data into wide format where each metric_name becomes a separate column
averages_wide <- averages %>%
  pivot_longer(cols = everything(), names_to = "full_colname", values_to = "average") %>%
  left_join(column_info, by = "full_colname") %>%
  mutate(participant_id = as.numeric(participant_id)) %>%
  select(participant_id, metric_name, average) %>%
  pivot_wider(names_from = metric_name, values_from = average)  # Pivot into wide format

# Step 5: Join with risk_data based on the participant_id value
final_data <- risk_data %>%
  left_join(averages_wide, by = c("idxP" = "participant_id"))
