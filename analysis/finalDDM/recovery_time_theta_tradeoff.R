# Do the latency (time) and non-decision time (theta) estimates trade off?
# Both are predicted by the same traits (Berlin +, GRiPS -) and correlate r = .20 across
# participants. If the model cannot fully separate them, part of that shared pattern could
# be an estimation artifact rather than two separate trait associations.
#
# Three checks, using the cached recovery fits (real scenario, b5000/s10000) and the
# original posterior; nothing is refit.
#   1. Within-participant posterior correlation of time[p] and theta.p[p] in the real fit.
#      Strongly negative = the data cannot tell a later onset from a longer non-decision time.
#   2. Correlation of recovery errors across participants. If errors in time and theta
#      move together, estimation couples the two parameters.
#   3. Between-participant correlation of time and theta, generating vs recovered. If
#      recovery inflates it, some of the observed r = .20 may be induced by the model.

main_folder <- "C:/Users/Brian Kim/Documents/OneDrive - Drexel University/Documents/classes/cognew/Neuroethics/IndividualDifferenceRiskBenefit/DATA_ANALYSIS/analysis/finalDDM"
out_dir <- file.path(main_folder, "Figures", "ParameterRecovery")
fit_suffix <- "_b5000_s10000"

load(file.path(main_folder, "FINAL_PARTICIPANT_LIST.RData"))
retained_idxP <- unique_idxP_table[[grep("idx", names(unique_idxP_table), ignore.case = TRUE)[1]]]

time_cols  <- sprintf("time[%d]", retained_idxP)
theta_cols <- sprintf("theta.p[%d]", retained_idxP)

# ---- 1. Within-participant posterior correlation (real data) ----
load(file.path(main_folder, "NeuroethicsPosteriorAndPredictions.RData"))   # m1_posterior
posterior_r <- mapply(function(tc, hc) cor(m1_posterior[[tc]], m1_posterior[[hc]]),
                      time_cols, theta_cols)
cat("1. Within-participant posterior r(time, theta), real fit\n")
print(round(quantile(posterior_r, c(.1, .25, .5, .75, .9)), 2))
cat("   proportion with r < -.5:", round(mean(posterior_r < -.5), 3), "\n\n")

# ---- 2 & 3. Recovery (real scenario) ----
load(file.path(out_dir, "simulated_data_real.RData"))                      # true_params
load(file.path(out_dir, sprintf("recovery_fit_real%s.RData", fit_suffix)))  # recovered_draws

true_time   <- true_params$time[retained_idxP]
true_theta  <- true_params$theta[retained_idxP]
recov_time  <- colMeans(recovered_draws[, time_cols])
recov_theta <- colMeans(recovered_draws[, theta_cols])

error_r <- cor(recov_time - true_time, recov_theta - true_theta)
cat("2. Correlation of recovery errors (recovered - true), time vs theta:", round(error_r, 3), "\n\n")

cat("3. Between-participant r(time, theta)\n")
cat("   generating values:", round(cor(true_time, true_theta), 3), "\n")
cat("   recovered values: ", round(cor(recov_time, recov_theta), 3), "\n")

write.csv(data.frame(idxP = retained_idxP, posterior_r, true_time, true_theta, recov_time, recov_theta),
          file.path(out_dir, sprintf("recovery_time_theta_tradeoff%s.csv", fit_suffix)), row.names = FALSE)
