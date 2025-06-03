# plotting_res_item_feature_importance.R - plotting, regression analysis for parameters
# to the items
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
# 2024/11/30      Kianté  Fernandez                   coded up version one

# Required libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)

# Load data
# load("/Users/kiante/Downloads/RidgeResults.RData")

category_colors <- c(
  "BIG5" = "#808080",
  "Moral Foundations (MFQ)" = "#009E73",
  "Religious & Spiritual" = "#F0E442", 
  "Risk Taking (GRiPS)" = "#d62728",
  "Interpersonal Reactivity" = "#88CCEE",
  "Values (Schwartz)" = "#AA4499",
  "Political" = "#117733",
  "Demographics" = "#882255",
  "Moral Commitments" = "#bcbd22",
  "Numeracy" = "#661100",
  "Optimize" = "#332288",
  "Human Mind" = "#D55E00"
  # "Other" = "#7f7f7f"
)

# Update category mapping to separate BIG5 dimensions based on scoring key
category_mapping <- list(
  "BIG5" = "^BIG5",
  
  # Rest of categories remain the same
  "Moral Foundations (MFQ)" = "^MFSQ",
  "Religious & Spiritual" = c("^Religion", "^AttendChurch", "Spirituality", "Religiosity", "^ReligiousActivity", "meditation"),
  "Risk Taking (GRiPS)" = "^GRiPS_",
  "Interpersonal Reactivity" = c("^PT", "^EC", "^PD", "^FS"),
  "Values (Schwartz)" = c("TRADITION", "SELF.DIRECTION", "SECURITY", "CONFORMITY", "BENEVOLENCE",
                          "POWER", "ACHIEVEMENT", "STIMULATION", "HEDONISM", "UNIVERSALISM"),
  "Political" = c("^Political_Party", "liberal", "conservative", "republican", "democrat", "independent"),
  "Demographics" = c("^education", "^RelationshipStatus", "age", "^sex", "race_ethnicity", "YearOfEducation"),
  "Moral Commitments" = "^MC_",
  "Numeracy" = c("^Berlin", "^SNS"),
  "Optimize" = "Optimize",
  "Human Mind" = "HumanMind"
)

# Define parameter names with expressions for rendering
measure_names <- c(
  "time" = expression(s[time]),
  "b_gain" = expression(nu[beta[gain]]),
  "b_risk" = expression(nu[beta[risk]]),
  "alpha" = expression(alpha),
  "theta" = expression(tau),
  "bias" = expression(z)
)

measure_names <- c(
  "time" = "Time",
  "b_gain" = "Gain Sensitivity",
  "b_risk" = "Risk Sensitivity",
  "alpha" = "Response Caution",
  "theta" = "Non-decision Time",
  "bias" = "Bias"
)

# Function to map to the categories we hard-coded above
get_category <- function(feature_name) {
  for(cat_name in names(category_mapping)) {
    patterns <- category_mapping[[cat_name]]
    if(any(sapply(patterns, function(p) grepl(p, feature_name)))) {
      return(cat_name)
    }
  }
  return("Other")
}

create_importance_plot <- function(data, measure) {
  
  data_plt <- data %>%
    mutate(
      Category = sapply(Feature, get_category),
      StdError = StdError,
      Importance = Importance
    ) %>%
    arrange(desc(Importance))
  
  # Get top 10 data for inset
  top_10_data <- data_plt %>%
    slice(1:10)
  
  data_plt$Category <- factor(data_plt$Category, levels = names(category_colors))
  # browser()
  # main plot
  main_plot <- ggplot(data_plt, 
                      aes(x = Importance, 
                          y = reorder(Feature, Importance),
                          fill = Category)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = category_colors, drop = FALSE) +
    scale_x_continuous(
      labels = scales::scientific_format(digits = 2),
      breaks = scales::pretty_breaks(n = 6),
      expand = c(0, 0)
    ) +
    theme_classic() +
    theme(
      text = element_text(family = "Helvetica", size = 18),           # to increased size base size
      axis.text = element_text(size = 16),                           # to increased size axis text
      axis.title = element_text(size = 20, face = "bold"),           # to increased size axis title
      axis.text.y = element_blank(),
      plot.title = element_text(size = 24, face = "bold",            # to increased size title size
                                margin = margin(b = 20)),
      legend.title = element_text(size = 20, face = "bold"),         # to increased size legend title
      legend.text = element_text(size = 16),                         # to increased size legend text
      legend.position = "right",
      legend.key.size = unit(1.5, "cm"),                            # to increased size legend key size
      axis.line.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")
    ) +
    labs(
      x = "Feature Importance",
      y = "items",
      fill = "traits",
      title = measure_names[measure]
    )
  
  # inset plot
  main_plot + 
    annotation_custom(
      ggplotGrob(
        ggplot(top_10_data, 
               aes(x = Importance, 
                   y = reorder(Feature, Importance),
                   fill = Category)) +
          geom_bar(stat = "identity", width = 1) +
          scale_fill_manual(values = category_colors, drop = FALSE) +
          scale_x_continuous(
            labels = scales::scientific_format(digits = 2)
          ) +
          theme_classic() +
          theme(
            text = element_text(family = "Helvetica", size = 16),          # to increased size
            axis.text.y = element_text(size = 14),                         # to increased size y-axis text
            axis.text.x = element_text(size = 14),                         # to increased size x-axis text
            legend.position = "none",
            panel.grid = element_blank(),
            plot.title = element_blank(),
            axis.title = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1.5)  # to increased size border
          )
      ),
      xmin = quantile(data_plt$Importance, 0.4),
      xmax = quantile(data_plt$Importance, 0.99),
      ymin = 10,
      ymax = 85
    )
}

# Create plots for each measure
measures <- c("time", "b_gain", "b_risk", "alpha", "theta", "bias")
plots <- lapply(measures, function(measure) {
  create_importance_plot(importance_results[[measure]], measure)
})

# Combine plots with a single legend
combined_plot <- wrap_plots(plots, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    plot.tag = element_text(size = 24, face = "bold")  # Control tag appearance here
  )

# Save combined plot
ggsave(
  filename = "Figures/feature_importance_combined.png",
  plot = combined_plot,
  width = 20,
  height = 24,
  dpi = 300
)

# Display combined plot
print(combined_plot)
