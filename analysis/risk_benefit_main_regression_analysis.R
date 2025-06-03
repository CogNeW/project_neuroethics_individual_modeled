# Fernandez_main_regression_analysis.R -

# Copyright (C) 2022 Kianté Fernandez, <kiantefernan@gmail.com>
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
# 2022/12/05      Kianté  Fernandez                     recoded analysis for communication

# Libraries

library(GGally) # Extension to 'ggplot2'
library(reshape2) # Flexibly Reshape Data: A Reboot of the Reshape Package
library(psych) # Procedures for Psychological, Psychometric, and Personality Research
library(car) # Companion to Applied Regression
library(knitr) # A General-Purpose Package for Dynamic Report Generation in R
library(broom) # Convert Statistical Objects into Tidy Tibbles
library(tidyverse) # Easily Install and Load the 'Tidyverse'
library(kableExtra) # Construct Complex Table with 'kable' and Pipe Syntax
library(scales) # Scale Functions for Visualization
library(outliers) # Tests for outliers
library(sjPlot) # Data Visualization for Statistics in Social Science
library(report) # Automated Reporting of Results and Statistical Models
library(emmeans) # Estimated Marginal Means, aka Least-Squares Means
library(lme4) # Linear Mixed-Effects Models using 'Eigen' and S4
library(lmerTest) # Tests in Linear Mixed Effects Models
library(pairwiseComparisons) # Multiple Pairwise Comparison Tests
library(here) # A Simpler Way to Find Your Files
library(doParallel) # Foreach Parallel Adaptor for the 'parallel' Package
library(ggpubr) # 'ggplot2' Based Publication Ready Plots
library(gghalves) # Compose Half-Half Plots Using Your Favourite Geoms

# gets last trial which has total task time and other covariates
source(here::here("ext", "get_last.R"))

# load clean data
Neuroethics_Judgement <- read.csv(here::here("data", "Cleaned_Pilot", "clean_data_longV2.csv"), stringsAsFactors = TRUE)
# select data used for analysis
data_h <- Neuroethics_Judgement %>%
  dplyr::select(user_id, condition, rt)
  # mutate(dominate = ifelse(Yes >= 75 * .90 | No >= 75 * .90, "dominate", "not_dominate"))
  # filter(!dominate == "dominate") %>%


## demographics
Neuroethics_Judgement %>%
  dplyr::select(user_id, trial_index, sex, age, race_ethnicity, education, Income, time2completion, time_elapsed, Yes, No) %>%
  mutate(
    time_elapsed = time_elapsed / 60000,
    trial_index = 1
  ) %>%
  group_by(user_id) %>%
  mutate(trial_index = cumsum(trial_index)) %>%
  mutate(education = fct_recode(education,
    HighSchool = "high school degree",
    Undergraduate = "bachelor",
    Postgraduate = "masters",
    Postgraduate = "Doctoral/Professional"
  )) %>%
  mutate(dominate = ifelse(Yes >= 75 * .90 | No >= 75 * .90, "dominate", "not_dominate")) %>%
  filter(!dominate == "dominate") -> demographics

demographics <- get_last(demographics)

# sex
n.female <- nrow(demographics[which(demographics$sex == "Female"), ])
n.male <- nrow(demographics[which(demographics$sex == "Male"), ])
n.other <- nrow(demographics[which(demographics$sex == "Other"), ])
# age
minage <- min(demographics$age, na.rm = T)
maxage <- max(demographics$age, na.rm = T)
mage <- round(mean(demographics$age, na.rm = T), 2)
sdage <- round(sd(demographics$age, na.rm = T), 2)
# race ethnicity
n.white <- nrow(demographics[which(demographics$race_ethnicity == "White"), ])
n.black <- nrow(demographics[which(demographics$race_ethnicity == "Black"), ])
n.asian <- nrow(demographics[which(demographics$race_ethnicity == "Asian"), ])
n.hispanic <- nrow(demographics[which(demographics$race_ethnicity == "Hispanic or Latino"), ])
n.othereth <- nrow(demographics[which(demographics$race_ethnicity == "Other"), ]) +
  nrow(demographics[which(demographics$race_ethnicity == "Mixed"), ])
n.pacific <- nrow(demographics[which(demographics$race_ethnicity == "Native Hawaiian/Pacific Islander"), ])
n.aindian <- nrow(demographics[which(demographics$race_ethnicity == "American Indian/Alaska Native"), ])
n.notreported <- nrow(demographics[which(demographics$race_ethnicity == "Prefer not to say"), ])
# education
n.HighSchool <- nrow(demographics[which(demographics$education == "HighSchool"), ])
n.Undergraduate <- nrow(demographics[which(demographics$education == "Undergraduate"), ])
n.Postgraduate <- nrow(demographics[which(demographics$education == "Postgraduate"), ])
# income
minincome <- min(demographics$Income, na.rm = T)
maxincome <- max(demographics$Income, na.rm = T)
medincome <- round(median(demographics$Income, na.rm = T), 2)
IQRincome <- round(IQR(demographics$Income, na.rm = T), 2)

# task
mtask <- round(mean(demographics$time_elapsed, na.rm = T), 2)
sdtask <- round(sd(demographics$time_elapsed, na.rm = T), 2)

# survey
msurvey <- round(mean(demographics$time2completion, na.rm = T), 2)
sdsurvey <- round(sd(demographics$time2completion, na.rm = T), 2)

# trials
mtrial <- round(mean(Neuroethics_Judgement$rt, na.rm = T), 2)
sdtrial <- round(sd(Neuroethics_Judgement$rt, na.rm = T), 2)

# yes
myes <- round(mean(demographics$Yes, na.rm = T), 0)
sdyes <- round(sd(demographics$Yes, na.rm = T), 0)

# No
mno <- round(mean(demographics$No, na.rm = T), 0)
sdno <- round(sd(demographics$No, na.rm = T), 0)

sample <- length(demographics$user_id)
excluded <- 616 - sample

exclusion_rate <- round(excluded / 616 * 100, 1)

paste(report::report_participants(demographics, spell_n = TRUE))

# create table
demographics %>%
  ungroup() %>%
  filter(!education == "less than high school") %>%
  mutate(education = factor(education)) %>%
  na.omit() %>%
  dplyr::select(sex, age, race_ethnicity, education, Income) %>%
  gtsummary::tbl_summary(by = education) %>%
  gtsummary::add_n() %>%
  gtsummary::as_gt()


GLMM_data <- Neuroethics_Judgement %>%
  select(user_id, condition, risk, gain, experimental_treatment_selected, rt) %>%
  mutate(experimental_treatment_selected = ifelse(experimental_treatment_selected == "Yes to treatment", 1, 0)) %>%
  group_by(user_id) %>%
  mutate(
    risk = scale(risk),
    gain = scale(gain)
  ) %>%
  ungroup()

# within subject scaling used here
# here is the first model
glmm.1 <- glmer(experimental_treatment_selected ~ risk + gain + (risk + gain | user_id),
data = GLMM_data, family = binomial(link = "logit"),
control = glmerControl(optimizer = "bobyqa"))
# extract the eblups for inference
eblups <- do.call(rbind.data.frame, ranef(glmm.1)) %>%
  rownames_to_column() %>%
  mutate(user_id = factor(str_remove(rowname, "user_id."))) %>%
  mutate(
    Intercept = `(Intercept)` + fixef(glmm.1)[["(Intercept)"]],
    risk = risk + fixef(glmm.1)[["risk"]],
    gain = gain + fixef(glmm.1)[["gain"]]
  ) %>%
  select(-rowname, -`(Intercept)`) %>%
  left_join(dplyr::select(data_h, user_id, condition), by = "user_id") %>%
  distinct()

se <- sqrt(diag(vcov(glmm.1)))
# table of estimates with 95% CI
(tab <- cbind(Est = fixef(glmm.1), LL = fixef(glmm.1) - 1.96 * se, UL = fixef(glmm.1) + 1.96 * se))
exp(tab)

# test for difference between parameters.
# (0.3134202) -  (-0.2497197)
car::linearHypothesis(glmm.1, "gain + risk = 0")

# test for differences between risk and gains
test <- eblups %>%
  select(user_id, condition, risk, gain) %>%
  reshape2::melt(id = c("user_id", "condition"), variable.name = "type", value.name = "subject_specific_effect") %>%
  mutate(
    user_id = factor(user_id),
    subject_specific_effect = abs(subject_specific_effect)
  )

t <- t.test(subject_specific_effect ~ type, paired = TRUE, alternative = "two.sided", conf.level = 0.95, var.equal = F, data = test)
tab <- map_df(list(t), tidy)
knitr::kable(tab[c(
  "estimate", "statistic", "p.value", "parameter", "conf.low",
  "conf.high", "method"
)], caption = "Parired T-test") %>%
  kable_styling(latex_options = "hold_position")

# plot
tt <- ggstatsplot::ggwithinstats(
  data = test,
  x = type,
  y = subject_specific_effect,
  title = "Paired t-test of the participant specific risk and benifit sensitivity",
  ylab = "Individual-level random slope coefficients",
  xlab = "sensitivity",
  subtitle = "",
  caption = "",
  ggstatsplot.layer = T,
  point.path = F,
  point.args = list(size = 3, alpha = 0.4, position = "jitter", width = 0.2),
  mean.plotting = T,
  effsize.type = "d",
  nboot = 1000L,
  bf.message = F,
)
tt
test %>%
  group_by(type) %>%
  summarise(
    median = median(subject_specific_effect),
    mean = mean(subject_specific_effect),
    sd = sd(subject_specific_effect)
  )

# test for interdoimain stability
test <- eblups %>%
  select(user_id, condition, risk, gain)

test %>%
  group_by(condition) %>%
  summarise(
    eBLUPs_risk = mean(risk),
    sd_risk = sd(risk),
    eBLUPs_gain = mean(gain),
    sd_gain = sd(gain)
  ) %>%
  kable(digits = 3) %>%
  kable_styling(latex_options = "hold_position")

aov.1 <- stats::oneway.test(risk ~ condition, data = test, var.equal = T)
aov.2 <- stats::oneway.test(gain ~ condition, data = test, var.equal = T)

aov_p1 <- ggstatsplot::ggbetweenstats(
  data = test,
  x = condition,
  y = risk,
  title = "Risk sensitivity across cognitive domains",
  ylab = "Individual-level random slope coefficients",
  xlab = "Cognitve Domain",
  pairwise.comparisons = F,
  bf.message = F
) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

aov_p2 <- ggstatsplot::ggbetweenstats(
  data = test,
  x = condition,
  y = gain,
  title = "Benefit sensitivity across cognitive domains",
  ylab = "Individual-level random slope coefficients",
  xlab = "Cognitve Domain",
  pairwise.comparisons = F,
  bf.message = F
) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# create a plot
aov_p1 + aov_p2

# create a table
tab <- map_df(list(aov.1, aov.2), tidy)
knitr::kable(tab[c("num.df", "den.df", "statistic", "p.value", "method")], caption = "one way anova testing cognitve domain; (Risk and Gain)") %>%
  kable_styling(latex_options = "hold_position")

# relevel so we compare to motor function
GLMM_data$condition <- relevel(GLMM_data$condition, "Motor Function")

glmm.cog <- glmer(experimental_treatment_selected ~ (risk + gain)*condition + (risk + gain | user_id),
data = GLMM_data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
summary(glmm.cog)

lmm.cog <- lmer(log(rt) ~ (risk + gain)*condition + (risk + gain | user_id),
                  data = GLMM_data, control = lmerControl(optimizer = "bobyqa"))
summary(lmm.cog)
