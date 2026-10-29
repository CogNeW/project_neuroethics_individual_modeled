# Supplemental Table 1: group-level recovery (real + centered scenarios, b5000/s10000 fits).
# Generating value = mean of the individual generating values (as in parameter_recovery.R Stage 4).

main_folder <- "C:/Users/Brian Kim/Documents/OneDrive - Drexel University/Documents/classes/cognew/Neuroethics/IndividualDifferenceRiskBenefit/DATA_ANALYSIS/analysis/finalDDM"
out_dir <- file.path(main_folder, "Figures", "ParameterRecovery")
fit_suffix <- "_b5000_s10000"

load(file.path(main_folder, "FINAL_PARTICIPANT_LIST.RData"))
retained_idxP <- unique_idxP_table[[grep("idx", names(unique_idxP_table), ignore.case = TRUE)[1]]]

params <- c(alpha = "Response caution", b_gain = "Cure-probability sensitivity",
            b_risk = "Harm-probability sensitivity", theta = "Non-decision time",
            time = "time", bias = "bias")
scenarios <- c(real = "Real", centered = "Centered")

rows <- list()
for (s in names(scenarios)) {
  load(file.path(out_dir, sprintf("simulated_data_%s.RData", s)))
  load(file.path(out_dir, sprintf("recovery_fit_%s%s.RData", s, fit_suffix)))
  for (p in names(params)) {
    mu <- recovered_draws[[paste0(p, ".mu")]]
    ci <- quantile(mu, c(.025, .975))
    gen <- mean(true_params[[p]][retained_idxP])
    rows[[length(rows) + 1]] <- data.frame(
      Scenario = scenarios[[s]], Parameter = params[[p]],
      Generating = gen, Recovered = mean(mu), CrI_lower = ci[[1]], CrI_upper = ci[[2]],
      Covered = gen >= ci[[1]] & gen <= ci[[2]], Rhat = rhat[paste0(p, ".mu"), 1])
  }
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(out_dir, sprintf("recovery_group_table%s.csv", fit_suffix)), row.names = FALSE)
print(format(tab, digits = 3), row.names = FALSE)
