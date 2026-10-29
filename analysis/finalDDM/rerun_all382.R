# ---------------------------------------------------------------------------
# Re-run the reported individual-difference analyses on all 382 fitted
# participants, instead of the 376 in FINAL_PARTICIPANT_LIST.RData
# ---------------------------------------------------------------------------
# Six participants (4BSARJqx, 5pd8VS5F, cLZrBIi5, pn5EKENt, R8pRzMWN, xq1MWUid)
# were fitted by the stDDM but dropped when unique_idxP_table was built. The code
# that built that list is not in this repository and no property of these six
# distinguishes them from the retained 376 -- see exclusion_robustness.R, which
# rules out survey completeness, attention checks, trial counts, choice profile,
# RT, repeat participation, shared IPs, geography, and an exhaustive search over
# 222 per-participant variables including every stDDM posterior mean and SD.
#
# Rather than justify an exclusion whose rule cannot be recovered, this script
# re-runs the analyses with those six added back and reports whether anything
# changes. The pipeline is copied from DDM_STDDM_Regression.Rmd; the ONLY
# difference is the participant filter, which is why that filter is a parameter
# here. Everything downstream -- the collapse of raw.csv, the factor coercions,
# the Berlin recode, the predictor screen, complete.cases, cv.glmnet lambda
# selection, and compute_feature_importance -- is reproduced verbatim.
#
# Note the published analyses do not actually run on 376: complete.cases on the
# predictor set drops a further two, giving 374 (which the manuscript reports).
# The comparison below therefore contrasts the published 374 against whatever
# the 382 set yields after the same complete.cases screen.
# ---------------------------------------------------------------------------

suppressMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(data.table); library(glmnet)
})

base_dir <- "C:/Users/Brian Kim/Documents/OneDrive - Drexel University/Documents/classes/cognew/Neuroethics/IndividualDifferenceRiskBenefit/DATA_ANALYSIS"
out_dir  <- file.path(base_dir, "analysis/finalDDM/Figures/ExclusionRobustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_lines <- character(0)
say <- function(...) { txt <- paste0(...); cat(txt, "\n", sep = ""); log_lines <<- c(log_lines, txt) }
show_obj <- function(x) { log_lines <<- c(log_lines, capture.output(print(x))); print(x) }
hdr <- function(x) say("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78))

# ---------------------------------------------------------------------------
# 1. Inputs
# ---------------------------------------------------------------------------
say("Loading posterior (large file, please wait)...")
ev <- new.env()
load(file.path(base_dir, "analysis/finalDDM/NeuroethicsPosteriorAndPredictions.RData"), envir = ev)
m1_posterior <- ev$m1_posterior
rm(ev); gc(verbose = FALSE)

e <- new.env(); load(file.path(base_dir, "data/FINAL_RT_EXCLUDED_DATA_20241014.RData"), envir = e)
risk_data <- e$risk_data %>% mutate(idxP = as.numeric(as.factor(user_id)))

el <- new.env(); load(file.path(base_dir, "analysis/finalDDM/FINAL_PARTICIPANT_LIST.RData"), envir = el)
ids_376 <- el$unique_idxP_table$user_id
ids_382 <- unique(risk_data$user_id)
the_six <- setdiff(ids_382, ids_376)

raw_data_all <- read.csv(file.path(base_dir, "data/raw_mturk_data/raw.csv"), stringsAsFactors = FALSE)

say("participants fitted by the stDDM: ", length(ids_382))
say("participants in FINAL_PARTICIPANT_LIST: ", length(ids_376))
say("added back here: ", paste(sort(the_six), collapse = ", "))

# ---------------------------------------------------------------------------
# 2. Per-participant stDDM parameters (verbatim from the Rmd)
# ---------------------------------------------------------------------------
selected_columns <- m1_posterior %>% select(matches("\\[\\d+\\]$"))
averages <- selected_columns %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
column_info <- tibble(
  full_colname   = colnames(averages),
  participant_id = str_extract(colnames(averages), "(?<=\\[)\\d+(?=\\])"),
  metric_name    = str_replace(colnames(averages), "\\.?p?\\[\\d+\\]$", "")
)
averages_wide <- averages %>%
  pivot_longer(cols = everything(), names_to = "full_colname", values_to = "average") %>%
  left_join(column_info, by = "full_colname") %>%
  mutate(participant_id = as.numeric(participant_id)) %>%
  select(participant_id, metric_name, average) %>%
  pivot_wider(names_from = metric_name, values_from = average)
rm(m1_posterior, selected_columns, averages); gc(verbose = FALSE)

risk_data_simplified <- risk_data %>% group_by(user_id) %>% slice(1) %>% ungroup()
averages_wide <- averages_wide %>%
  left_join(risk_data_simplified %>% select(idxP, user_id), by = c("participant_id" = "idxP"))

# ---------------------------------------------------------------------------
# 3. Build final_table for an arbitrary participant set (verbatim from the Rmd)
# ---------------------------------------------------------------------------
build_final_table <- function(keep_ids) {
  raw_data <- raw_data_all %>% filter(user_id %in% keep_ids)   # <- the only difference
  collapsed_data <- raw_data %>% group_by(user_id) %>% summarise(across(everything(), ~ first(.)))

  ft <- left_join(collapsed_data, averages_wide, by = "user_id")
  ft <- as.data.table(ft)

  columns_to_factor <- c(
    "sex", "education", "ReligionRaised", "ReligionNow", "AttendChurch",
    "ReligiousActivity", "meditation", "RelationshipStatus",
    "Major", "EmploymentStatus", "Political_Party", "Optimize", "HumanMind"
  )
  risk_data_subset <- risk_data_simplified[, c("user_id", columns_to_factor)]
  ft <- ft %>% select(-all_of(columns_to_factor))
  ft <- merge(ft, risk_data_subset, by = "user_id")

  ft[YearOfEducation == "two", YearOfEducation := "2"]
  ft[, YearOfEducation := gsub("[^0-9.]", "", YearOfEducation)]
  ft[, YearOfEducation := as.numeric(YearOfEducation)]

  ft$Berlin1 <- factor(ifelse(ft$Berlin1 == 25, 1, 0))
  ft$Berlin2 <- factor(ifelse(ft$Berlin2 == 30, 1, 0))
  ft$Berlin3 <- factor(ifelse(ft$Berlin3 == 20, 1, 0))
  ft$Berlin4 <- factor(ifelse(ft$Berlin4 == 50, 1, 0))
  for (v in c("MC_Lying", "MC_Assassin", "MC_Torture", "MC_Murder", "MC_Stealing", "MC_Steril")) {
    ft[[v]] <- factor(ft[[v]])
  }
  ft
}

# compute_feature_importance, verbatim from DDM_STDDM_Regression.Rmd
compute_feature_importance <- function(data, response, predictors,
        n_subsets = 100, subset_fraction = 1, feature_fraction = 0.8, lambda = 0.1) {
  independent_formula <- as.formula(paste("~", paste(predictors, collapse = " + ")))
  full_matrix <- model.matrix(independent_formula, data = data)[, -1]
  all_columns <- colnames(full_matrix)
  feature_frequencies <- colMeans(full_matrix != 0)
  feature_importance <- matrix(NA, nrow = n_subsets, ncol = length(all_columns))
  colnames(feature_importance) <- all_columns
  for (i in 1:n_subsets) {
    subset_data <- data %>% sample_frac(subset_fraction)
    data_matrix <- model.matrix(independent_formula, data = subset_data)[, -1]
    selected_features <- sample(all_columns, size = floor(feature_fraction * length(all_columns)))
    data_matrix <- data_matrix[, selected_features, drop = FALSE]
    missing_cols <- setdiff(all_columns, colnames(data_matrix))
    for (col in missing_cols) {
      data_matrix <- cbind(data_matrix, 0); colnames(data_matrix)[ncol(data_matrix)] <- col
    }
    data_matrix <- data_matrix[, all_columns, drop = FALSE]
    y <- subset_data[[response]]
    ridge_model <- glmnet(data_matrix, y, alpha = 0, lambda = lambda, standardize = TRUE)
    feature_importance[i, ] <- abs(as.vector(coef(ridge_model)[-1]))
  }
  importance_scores <- colMeans(feature_importance, na.rm = TRUE)
  importance_se <- apply(feature_importance, 2, function(x) sd(x, na.rm = TRUE) / sqrt(n_subsets))
  data.frame(Feature = colnames(feature_importance),
             Importance = importance_scores * feature_frequencies,
             StdError = importance_se) %>% arrange(desc(Importance))
}

dependent_vars <- c("time", "b_gain", "b_risk", "alpha", "theta", "bias")

run_pipeline <- function(keep_ids, label, seed = 1234) {
  set.seed(seed)   # 1234 as in the Rmd; varied to measure the method's own noise
  ft <- build_final_table(keep_ids)
  column_names <- colnames(ft)
  predictors <- column_names[
    !grepl("^(participant_id|user_id|rt|key_press|trial_index|time_elapsed|condition|condition_per_trial|total_trials|treatment_display|risk|gain|no_effect|experimental_treatment_selected|Progress|Duration..in.seconds.|RecordedDate|ResponseId|attention_check|Q1802|Income|Major|EmploymentStatus)", column_names) &
    column_names != "X" &
    !column_names %in% c("MFSQ6", "MFSQ22")]   # MFQ catch items, as in the main Rmd
  ft <- ft[complete.cases(ft[, ..predictors])]
  predictors <- setdiff(predictors, dependent_vars)

  independent_formula <- as.formula(paste("~", paste(predictors, collapse = " + ")))
  data_matrix <- model.matrix(independent_formula, data = ft)[, -1]

  res <- list(n = nrow(ft), ids = ft$user_id, importance = list(), fit = data.frame())
  for (dv in dependent_vars) {
    y <- scale(ft[[dv]])
    cv_fit <- cv.glmnet(data_matrix, y, alpha = 0)
    best_lambda <- cv_fit$lambda.min
    rm_ <- glmnet(data_matrix, y, alpha = 0, lambda = best_lambda)
    pred <- predict(rm_, newx = data_matrix)
    r2 <- 1 - sum((y - pred)^2) / sum((y - mean(y))^2)
    res$fit <- rbind(res$fit, data.frame(dv = dv, lambda = best_lambda, r2_pct = 100 * r2))
    res$importance[[dv]] <- compute_feature_importance(ft, dv, predictors,
                                                       n_subsets = 100, lambda = best_lambda)
  }
  say("[", label, "] n after complete.cases = ", res$n)
  res
}

# ---------------------------------------------------------------------------
# 4. Run both
# ---------------------------------------------------------------------------
hdr("Re-running the ridge feature-importance analysis")
pub <- run_pipeline(ids_376, "published 376 list")
all <- run_pipeline(ids_382, "all 382 fitted")

say("\nParticipants actually analysed:")
say("  published pipeline: ", pub$n, "   (manuscript reports 374)")
say("  all-382 pipeline:   ", all$n)
say("  difference:         ", all$n - pub$n)

# CONTROL: the same 376 participants, only a different random seed. Both
# cv.glmnet and compute_feature_importance are stochastic (100 random subsets,
# each using a random 80% of features), so the importance ranking moves from run
# to run even with the data held fixed. Without this control, any movement
# between the 376 and 382 results would be over-read as an effect of the six.
ctl <- run_pipeline(ids_376, "published 376 list, different seed", seed = 99)

hdr("Model fit per dependent variable")
cmp_fit <- merge(pub$fit, all$fit, by = "dv", suffixes = c("_376", "_382"))
cmp_fit$d_r2 <- round(cmp_fit$r2_pct_382 - cmp_fit$r2_pct_376, 3)
show_obj(cmp_fit)

hdr("Do the feature-importance results change?")
summ <- data.frame()
for (dv in dependent_vars) {
  a <- pub$importance[[dv]]; b <- all$importance[[dv]]
  m <- merge(a, b, by = "Feature", suffixes = c("_376", "_382"))
  summ <- rbind(summ, data.frame(
    dv = dv,
    n_features = nrow(m),
    pearson_r  = round(cor(m$Importance_376, m$Importance_382), 4),
    spearman_r = round(cor(m$Importance_376, m$Importance_382, method = "spearman"), 4),
    top10_overlap = length(intersect(head(a$Feature, 10), head(b$Feature, 10))),
    top1_376 = a$Feature[1], top1_382 = b$Feature[1]))
}
show_obj(summ)

hdr("Top 10 features per dependent variable, side by side")
for (dv in dependent_vars) {
  say("\n--- ", dv, " ---")
  side <- data.frame(
    rank = 1:10,
    published_376 = head(pub$importance[[dv]]$Feature, 10),
    all_382       = head(all$importance[[dv]]$Feature, 10))
  show_obj(side)
}

hdr("CONTROL: same 376 people, different seed")
say("If adding the six moves the results no more than simply re-running the")
say("identical analysis with a different random seed, then the six make no")
say("detectable difference and the exclusion is immaterial.")
ctl_summ <- data.frame()
for (dv in dependent_vars) {
  a <- pub$importance[[dv]]; b <- ctl$importance[[dv]]; c2 <- all$importance[[dv]]
  m1 <- merge(a, b, by = "Feature", suffixes = c("_a", "_b"))
  m2 <- merge(a, c2, by = "Feature", suffixes = c("_a", "_b"))
  ctl_summ <- rbind(ctl_summ, data.frame(
    dv = dv,
    r_376_vs_376_newseed = round(cor(m1$Importance_a, m1$Importance_b, method = "spearman"), 4),
    r_376_vs_382         = round(cor(m2$Importance_a, m2$Importance_b, method = "spearman"), 4),
    top10_overlap_seed   = length(intersect(head(a$Feature, 10), head(b$Feature, 10))),
    top10_overlap_382    = length(intersect(head(a$Feature, 10), head(c2$Feature, 10)))))
}
show_obj(ctl_summ)
say("\nMean Spearman r, 376 vs the same 376 under a new seed: ",
    sprintf("%.4f", mean(ctl_summ$r_376_vs_376_newseed)))
say("Mean Spearman r, 376 vs all 382:                        ",
    sprintf("%.4f", mean(ctl_summ$r_376_vs_382)))
say("Mean top-10 overlap, new seed: ", sprintf("%.1f", mean(ctl_summ$top10_overlap_seed)), " of 10")
say("Mean top-10 overlap, all 382:  ", sprintf("%.1f", mean(ctl_summ$top10_overlap_382)), " of 10")

hdr("Rank movement of every feature (how far anything shifts)")
mv <- data.frame()
for (dv in dependent_vars) {
  a <- pub$importance[[dv]]; b <- all$importance[[dv]]; c2 <- ctl$importance[[dv]]
  a$r1 <- seq_len(nrow(a)); b$r2 <- seq_len(nrow(b)); c2$r3 <- seq_len(nrow(c2))
  m <- merge(merge(a[, c("Feature", "r1")], b[, c("Feature", "r2")], by = "Feature"),
             c2[, c("Feature", "r3")], by = "Feature")
  mv <- rbind(mv, data.frame(
    dv = dv,
    median_shift_382     = median(abs(m$r1 - m$r2)),
    median_shift_newseed = median(abs(m$r1 - m$r3)),
    max_shift_382        = max(abs(m$r1 - m$r2)),
    max_shift_newseed    = max(abs(m$r1 - m$r3))))
}
show_obj(mv)
say("\nRead the columns in pairs: rank movement from adding the six versus rank")
say("movement from nothing but a new random seed on identical data.")

for (dv in dependent_vars) {
  write.csv(pub$importance[[dv]], file.path(out_dir, paste0("importance_376_", dv, ".csv")), row.names = FALSE)
  write.csv(all$importance[[dv]], file.path(out_dir, paste0("importance_382_", dv, ".csv")), row.names = FALSE)
}

hdr("Conclusion")
say("The six participants were added back and the published pipeline re-run")
say("unchanged in every other respect, alongside a control that re-runs the")
say("published 376 under a different random seed.")
say("")
say("The ridge feature-importance ranking is itself stochastic, so the honest")
say("comparison is not 'did anything move' but 'did adding the six move it more")
say("than the method moves on its own'. The control answers that directly.")

writeLines(log_lines, file.path(out_dir, "rerun_all382_report.txt"))
say("\nWritten: ", file.path(out_dir, "rerun_all382_report.txt"))
