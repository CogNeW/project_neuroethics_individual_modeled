# Supplemental Figure 2: stDDM parameter recovery (real + centered scenarios, b5000/s10000 fits).
# Reads the cached outputs of parameter_recovery.R; does not refit anything.

library(dplyr)
library(ggplot2)

main_folder <- "C:/Users/Brian Kim/Documents/OneDrive - Drexel University/Documents/classes/cognew/Neuroethics/IndividualDifferenceRiskBenefit/DATA_ANALYSIS/analysis/finalDDM"
out_dir <- file.path(main_folder, "Figures", "ParameterRecovery")
fit_suffix <- "_b5000_s10000"

load(file.path(main_folder, "FINAL_PARTICIPANT_LIST.RData"))  # unique_idxP_table
retained_idxP <- unique_idxP_table[[grep("idx", names(unique_idxP_table), ignore.case = TRUE)[1]]]

param_prefix <- c(alpha = "alpha.p[%d]", b_gain = "b_gain.p[%d]", b_risk = "b_risk.p[%d]",
                  theta = "theta.p[%d]", time = "time[%d]", bias = "bias[%d]")
param_labels <- c(alpha = "Response caution", b_gain = "Cure-probability sensitivity",
                  b_risk = "Harm-probability sensitivity", theta = "Non-decision time",
                  time = "time", bias = "bias")
scenario_labels <- c(real = "Recovery", centered = "Centered")

load_scenario <- function(scenario) {
  load(file.path(out_dir, sprintf("simulated_data_%s.RData", scenario)))            # true_params
  load(file.path(out_dir, sprintf("recovery_fit_%s%s.RData", scenario, fit_suffix)))  # recovered_draws
  bind_rows(lapply(names(param_prefix), function(par) {
    d <- recovered_draws[, sprintf(param_prefix[[par]], retained_idxP)]
    data.frame(scenario = scenario, parameter = par, idxP = retained_idxP,
               true = true_params[[par]][retained_idxP], recovered = colMeans(d))
  }))
}

recovery_long <- bind_rows(lapply(names(scenario_labels), load_scenario))
write.csv(recovery_long, file.path(out_dir, sprintf("recovery_long%s.csv", fit_suffix)), row.names = FALSE)

plot_data <- recovery_long %>%
  mutate(parameter = factor(param_labels[parameter], levels = param_labels),
         scenario  = factor(scenario_labels[scenario], levels = scenario_labels))

# Same x and y range within each panel, so the identity line sits on the diagonal.
panel_range <- plot_data %>%
  group_by(scenario, parameter) %>%
  summarise(lo = min(true, recovered), hi = max(true, recovered), .groups = "drop") %>%
  tidyr::pivot_longer(c(lo, hi), values_to = "v")

fig <- ggplot(plot_data, aes(true, recovered)) +
  geom_blank(data = panel_range, aes(v, v)) +
  geom_abline(aes(slope = 1, intercept = 0, linetype = "Perfect recovery"),
              color = "grey30", linewidth = 0.7) +
  geom_point(aes(color = scenario), alpha = 0.45, size = 1.3) +
  # facet_wrap (not facet_grid) so every panel gets its own y-axis; color marks the row.
  facet_wrap(~ scenario + parameter, nrow = 2, scales = "free",
             labeller = function(labels) label_wrap_gen(16)(labels["parameter"])) +
  scale_color_manual(values = c(Recovery = "#2C7BB6", Centered = "#D7642C"), name = "Scenario") +
  scale_linetype_manual(values = c("Perfect recovery" = "dashed"), name = NULL) +
  guides(color = guide_legend(order = 1, override.aes = list(alpha = 1, size = 3)),
         linetype = guide_legend(order = 2)) +
  labs(x = "Generating value", y = "Recovered posterior mean",
       title = "stDDM parameter recovery") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        strip.text = element_text(face = "bold", size = 13),
        strip.background = element_rect(fill = "grey92", color = NA),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 13))

ggsave(file.path(out_dir, sprintf("recovery_scatter_combined%s.png", fit_suffix)),
       fig, width = 16, height = 7.5, dpi = 300)
