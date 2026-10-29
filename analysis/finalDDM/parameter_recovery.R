# Parameter recovery for the stDDM (M1) latency parameter
#
# Question: are individual latency (time[p]) estimates identifiable at our trial counts, or is
# the ~90% "harm-first" classification partly an artifact of shrinkage toward a positive
# group mean?
#
# Stages (each caches its output, so the script can be re-run and will skip finished stages):
#   1. Posterior widths of the real individual latency estimates (no fitting needed)
#   2. Simulate choice/RT data from the time-varying DDM using each participant's real trials
#      and known "true" parameters, under two scenarios:
#        "real"     - true params = each participant's posterior means (39/376 negative latency)
#        "centered" - same, but latency shifted so its mean is 0 (about half negative).
#                     If the model were just echoing a positive prior/group mean, recovered
#                     latencies here would still come out mostly positive.
#   3. Refit M1 in JAGS to each simulated data set with the original settings
#   4. Compare true vs recovered parameters; save tables and figures
#
# Usage (from DATA_ANALYSIS/analysis/finalDDM):
#   Rscript parameter_recovery.R            # runs every stage, both scenarios
#   Rscript parameter_recovery.R real       # one scenario only
# Stage 3 is the slow part (a full hierarchical fit per scenario); the two scenarios can be
# run in separate R sessions at the same time.
#
# SIGN CONVENTION: in stddm_model_M1.txt, time > 0 means the harm (risk) term enters first and
# the cure (gain) term is delayed. In 2ddm_r_cpp_M1.cpp, tIn > 0 delays the *harm* term
# (vec_tHealth). The simulator is therefore called with tIn = -time.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(runjags)
  library(coda)
  library(Rcpp)
})

main_folder <- "C:/Users/Brian Kim/Documents/OneDrive - Drexel University/Documents/classes/cognew/Neuroethics/IndividualDifferenceRiskBenefit/DATA_ANALYSIS/analysis/finalDDM"
behavior_path <- "C:/Users/Brian Kim/Documents/OneDrive - Drexel University/Documents/classes/cognew/Neuroethics/IndividualDifferenceRiskBenefit/DATA_ANALYSIS/data/FINAL_RT_EXCLUDED_DATA_20241014.RData"
out_dir <- file.path(main_folder, "Figures", "ParameterRecovery")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#   Rscript parameter_recovery.R real 5000 10000   # scenario, burnin, sample
args <- commandArgs(trailingOnly = TRUE)
scenarios <- if (length(args) >= 1) args[1] else c("real", "centered")
n_burnin  <- if (length(args) >= 2) as.integer(args[2]) else 1000
n_sample  <- if (length(args) >= 3) as.integer(args[3]) else 5000
# Fits with the original length keep their original file names; longer runs get a suffix so
# they don't overwrite (or get mistaken for) the earlier cached fits.
fit_suffix <- if (n_burnin == 1000 && n_sample == 5000) "" else sprintf("_b%d_s%d", n_burnin, n_sample)

load(file.path(main_folder, "NeuroethicsPosteriorAndPredictions.RData"))  # m1_posterior
load(file.path(main_folder, "FINAL_PARTICIPANT_LIST.RData"))              # unique_idxP_table
retained_idxP <- unique_idxP_table[[grep("idx", names(unique_idxP_table), ignore.case = TRUE)[1]]]

load(behavior_path)
risk_data <- risk_data %>%
  mutate(
    rt = rt / 1000,
    response = ifelse(experimental_treatment_selected == "Yes to treatment", 1, 0),
    risk = risk * 0.01,
    gain = gain * 0.01,
    idxP = as.numeric(as.factor(user_id))
  )
n_participants <- length(unique(risk_data$idxP))


# ---- Stage 1: posterior widths of the real latency estimates --------------------------------

time_cols <- paste0("time[", seq_len(n_participants), "]")
time_draws <- as.matrix(as.data.frame(m1_posterior)[, time_cols])

real_time_summary <- data.frame(
  idxP     = seq_len(n_participants),
  mean     = colMeans(time_draws),
  sd       = apply(time_draws, 2, sd),
  lower_95 = apply(time_draws, 2, quantile, 0.025),
  upper_95 = apply(time_draws, 2, quantile, 0.975),
  p_pos    = colMeans(time_draws > 0)
) %>%
  filter(idxP %in% retained_idxP) %>%
  mutate(
    classification = case_when(
      lower_95 > 0 ~ "Harm first (95% CrI > 0)",
      upper_95 < 0 ~ "Cure first (95% CrI < 0)",
      TRUE         ~ "Uncertain (CrI includes 0)"
    )
  )

cat("\n==== Stage 1: individual latency posteriors (n =", nrow(real_time_summary), ") ====\n")
cat("Posterior mean > 0:", sum(real_time_summary$mean > 0),
    "  < 0:", sum(real_time_summary$mean < 0), "\n")
print(table(real_time_summary$classification))
cat("Posterior SD: median", round(median(real_time_summary$sd), 3),
    "IQR", paste(round(quantile(real_time_summary$sd, c(.25, .75)), 3), collapse = "-"), "s\n")
cat("Group-level SD of posterior means:", round(sd(real_time_summary$mean), 3), "s\n")
cat("P(time > 0) >= .95:", sum(real_time_summary$p_pos >= .95),
    "  <= .05:", sum(real_time_summary$p_pos <= .05), "\n")

caterpillar <- real_time_summary %>%
  arrange(mean) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(x = rank, y = mean, ymin = lower_95, ymax = upper_95, color = classification)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = mean(m1_posterior[["time.mu"]]), color = "grey40") +
  geom_linerange(alpha = 0.6) +
  geom_point(size = 0.8) +
  scale_color_manual(values = c("Harm first (95% CrI > 0)" = "#D55E00",
                                "Cure first (95% CrI < 0)" = "#0072B2",
                                "Uncertain (CrI includes 0)" = "grey60")) +
  labs(x = "Participant (sorted by posterior mean)", y = "Latency (s)", color = NULL,
       title = "Individual latency posteriors (mean and 95% credible interval)",
       subtitle = "Grey line: group-level mean; dashed line: 0") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "latency_posterior_caterpillar.png"), caterpillar,
       width = 10, height = 6, dpi = 300)
write.csv(real_time_summary, file.path(out_dir, "latency_posterior_widths.csv"), row.names = FALSE)


# ---- Stage 2: simulate data --------------------------------------------------------------

sourceCpp(file.path(main_folder, "2ddm_r_cpp_M1.cpp"))

posterior_means <- function(prefix) {
  sapply(seq_len(n_participants), function(p) mean(m1_posterior[[sprintf(prefix, p)]]))
}
true_params_real <- data.frame(
  idxP   = seq_len(n_participants),
  b_gain = posterior_means("b_gain.p[%d]"),
  b_risk = posterior_means("b_risk.p[%d]"),
  alpha  = posterior_means("alpha.p[%d]"),
  theta  = posterior_means("theta.p[%d]"),
  time   = posterior_means("time[%d]"),
  bias   = posterior_means("bias[%d]")
)

make_true_params <- function(scenario) {
  tp <- true_params_real
  if (scenario == "centered") tp$time <- tp$time - mean(tp$time)
  tp
}

simulate_scenario <- function(scenario) {
  sim_file <- file.path(out_dir, sprintf("simulated_data_%s.RData", scenario))
  if (file.exists(sim_file)) {
    cat("Stage 2 (", scenario, "): using cached ", sim_file, "\n", sep = "")
    load(sim_file)
    return(list(sim_data = sim_data, true_params = true_params))
  }
  true_params <- make_true_params(scenario)
  sim_data <- risk_data %>%
    select(idxP, gain, risk) %>%
    left_join(true_params, by = "idxP")

  # One simulated trial per real trial (same displayed probabilities).
  # tIn = -time: see SIGN CONVENTION at the top.
  sim_data$y <- mapply(
    function(b_gain, b_risk, alpha, theta, time, bias, gain, risk) {
      ddm2_parallel(b_gain, b_risk, alpha, theta, -time, bias, gain, risk, 1, 1)
    },
    sim_data$b_gain, sim_data$b_risk, sim_data$alpha, sim_data$theta,
    sim_data$time, sim_data$bias, sim_data$gain, sim_data$risk
  )
  sim_data <- sim_data %>%
    mutate(rt = abs(y), response = as.integer(y > 0)) %>%
    select(idxP, gain, risk, y, rt, response)

  cat("Stage 2 (", scenario, "): simulated ", nrow(sim_data), " trials; ",
      "endorse rate ", round(mean(sim_data$response), 3),
      " (real ", round(mean(risk_data$response), 3), "); ",
      "median RT ", round(median(sim_data$rt), 2), " s (real ", round(median(risk_data$rt), 2), "); ",
      "trials hitting the 20 s cap: ", sum(sim_data$rt >= 20), "\n", sep = "")
  save(sim_data, true_params, file = sim_file)
  list(sim_data = sim_data, true_params = true_params)
}


# ---- Stage 3: refit M1 ---------------------------------------------------------------------

fit_scenario <- function(scenario, sim_data) {
  fit_file <- file.path(out_dir, sprintf("recovery_fit_%s%s.RData", scenario, fit_suffix))
  if (file.exists(fit_file)) {
    cat("Stage 3 (", scenario, "): using cached ", fit_file, "\n", sep = "")
    load(fit_file)
    return(recovered_draws)
  }

  # Same data format, inits, and MCMC settings as kim_analysis_stddm.R
  dat <- with(sim_data, dump.format(list(
    N = length(y), y = y, idxP = idxP, gain = gain, risk = risk, rt = rt,
    ns = length(unique(idxP))
  )))
  inits <- list(
    dump.format(list(alpha.mu = 2.4, alpha.pr = 0.05, theta.mu = 0.15, theta.pr = 0.05,
                     b_gain.mu = 0.1, b_gain.pr = 0.05, b_risk.mu = 0.05, b_risk.pr = 0.05,
                     bias.mu = 0.3, bias.kappa = 1,
                     .RNG.name = "base::Mersenne-Twister", .RNG.seed = 6666)),
    dump.format(list(alpha.mu = 2.2, alpha.pr = 0.05, theta.mu = 0.01, theta.pr = 0.05,
                     b_gain.mu = 0.3, b_gain.pr = 0.05, b_risk.mu = 0.1, b_risk.pr = 0.05,
                     bias.mu = 0.5, bias.kappa = 1,
                     .RNG.name = "base::Wichmann-Hill", .RNG.seed = 1234)),
    dump.format(list(alpha.mu = 2, alpha.pr = 0.05, theta.mu = 0.1, theta.pr = 0.05,
                     b_gain.mu = 0.3, b_gain.pr = 0.05, b_risk.mu = 0.01, b_risk.pr = 0.05,
                     bias.mu = 0.4, bias.kappa = 1,
                     .RNG.name = "base::Super-Duper", .RNG.seed = 99999))
  )
  monitor <- c("time.mu", "b_gain.mu", "b_risk.mu", "alpha.mu", "theta.mu", "bias.mu",
               "time", "b_gain.p", "b_risk.p", "alpha.p", "theta.p", "bias")

  cat("Stage 3 (", scenario, "): fitting M1 in JAGS, started ", format(Sys.time()), "\n", sep = "")
  fit <- run.jags(
    model = file.path(main_folder, "stddm_model_M1.txt"),
    monitor = monitor, data = dat, n.chains = 3, inits = inits,
    plots = FALSE, method = "parallel", summarise = FALSE, module = "wiener",
    burnin = n_burnin, sample = n_sample, thin = 1
  )
  samples <- as.mcmc.list(fit)
  recovered_draws <- as.data.frame(do.call(rbind, lapply(samples, as.matrix)))
  rhat <- gelman.diag(samples, multivariate = FALSE, autoburnin = FALSE)$psrf
  cat("Stage 3 (", scenario, "): finished ", format(Sys.time()),
      " (burnin ", n_burnin, ", sample ", n_sample, ")\n", sep = "")
  is_group <- endsWith(rownames(rhat), ".mu")
  print(round(rhat[is_group, ], 3))
  cat("Individual parameters:", sum(!is_group),
      "  R-hat > 1.01:", sum(rhat[!is_group, 1] > 1.01),
      "  > 1.05:", sum(rhat[!is_group, 1] > 1.05),
      "  max:", round(max(rhat[!is_group, 1]), 3), "\n")
  save(recovered_draws, rhat, file = fit_file)
  recovered_draws
}


# ---- Stage 4: compare true vs recovered ----------------------------------------------------

param_prefix <- c(time = "time[%d]", b_gain = "b_gain.p[%d]", b_risk = "b_risk.p[%d]",
                  alpha = "alpha.p[%d]", theta = "theta.p[%d]", bias = "bias[%d]")

compare_scenario <- function(scenario, true_params, recovered_draws) {
  long <- bind_rows(lapply(names(param_prefix), function(par) {
    cols <- sprintf(param_prefix[[par]], retained_idxP)
    d <- recovered_draws[, cols]
    data.frame(
      scenario  = scenario,
      parameter = par,
      idxP      = retained_idxP,
      true      = true_params[[par]][retained_idxP],
      recovered = colMeans(d),
      lower_95  = apply(d, 2, quantile, 0.025),
      upper_95  = apply(d, 2, quantile, 0.975)
    )
  })) %>%
    mutate(covered = true >= lower_95 & true <= upper_95)

  stats <- long %>%
    group_by(scenario, parameter) %>%
    summarise(
      r              = cor(true, recovered),
      shrinkage_slope = unname(coef(lm(recovered ~ true))[2]),
      rmse           = sqrt(mean((recovered - true)^2)),
      coverage_95    = mean(covered),
      .groups = "drop"
    )

  lat <- long %>% filter(parameter == "time")
  sign_table <- table(
    true      = ifelse(lat$true > 0, "harm first", "cure first"),
    recovered = ifelse(lat$recovered > 0, "harm first", "cure first")
  )
  cat("\n==== Stage 4 (", scenario, ") ====\n", sep = "")
  print(stats, n = Inf)
  cat("\nLatency sign: true vs recovered posterior mean\n")
  print(sign_table)
  cat("Sign agreement:", round(mean(sign(lat$true) == sign(lat$recovered)), 3),
      "  True harm-first:", sum(lat$true > 0), "/", nrow(lat),
      "  Recovered harm-first:", sum(lat$recovered > 0), "/", nrow(lat), "\n")
  cat("Recovered classification (95% CrI):",
      "harm first", sum(lat$lower_95 > 0),
      "| cure first", sum(lat$upper_95 < 0),
      "| uncertain", sum(lat$lower_95 <= 0 & lat$upper_95 >= 0), "\n")
  cat("Group-level time.mu: true mean of individual values", round(mean(true_params$time), 3),
      " recovered", round(mean(recovered_draws[["time.mu"]]), 3),
      " 95% CrI [", paste(round(quantile(recovered_draws[["time.mu"]], c(.025, .975)), 3),
                        collapse = ", "), "]\n")

  list(long = long, stats = stats, sign_table = sign_table)
}

results <- list()
for (scenario in scenarios) {
  sim <- simulate_scenario(scenario)
  draws <- fit_scenario(scenario, sim$sim_data)
  results[[scenario]] <- compare_scenario(scenario, sim$true_params, draws)
}

recovery_long  <- bind_rows(lapply(results, `[[`, "long"))
recovery_stats <- bind_rows(lapply(results, `[[`, "stats"))
write.csv(recovery_stats, file.path(out_dir, sprintf("recovery_stats_%s%s.csv", paste(scenarios, collapse = "_"), fit_suffix)),
          row.names = FALSE)

scatter <- ggplot(recovery_long, aes(true, recovered)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.4, size = 0.8) +
  facet_grid(scenario ~ parameter, scales = "free") +
  labs(x = "True value", y = "Recovered posterior mean",
       title = "stDDM parameter recovery (dashed line = perfect recovery)") +
  theme_minimal()
ggsave(file.path(out_dir, sprintf("recovery_scatter_%s%s.png", paste(scenarios, collapse = "_"), fit_suffix)),
       scatter, width = 14, height = 3 + 2.5 * length(scenarios), dpi = 300)

latency_plot <- recovery_long %>%
  filter(parameter == "time") %>%
  ggplot(aes(true, recovered, ymin = lower_95, ymax = upper_95)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_vline(xintercept = 0, color = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_linerange(alpha = 0.25) +
  geom_point(alpha = 0.6, size = 0.9) +
  facet_wrap(~ scenario) +
  labs(x = "True latency (s)", y = "Recovered latency (s), mean and 95% CrI",
       title = "Latency recovery",
       subtitle = "Off-diagonal quadrants = sign misclassified") +
  theme_minimal()
ggsave(file.path(out_dir, sprintf("latency_recovery_%s%s.png", paste(scenarios, collapse = "_"), fit_suffix)),
       latency_plot, width = 5 + 4 * length(scenarios), height = 5, dpi = 300)
