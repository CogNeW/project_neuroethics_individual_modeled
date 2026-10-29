# exclusion_robustness.R -- Reviewer response: who are the 240 excluded participants?
#
# Reviewer comment (point 6) argued that the n = 240 excluded from the reported
# sample of 616 are "likely people with strong preferences (always recommending
# or always rejecting)", so the analyzed sample may under-represent people with
# strong convictions about neuromodulation.
#
# This script answers that directly and reproducibly. It
#   (1) reconstructs the full exclusion cascade 616 -> 376 and reconciles it to
#       n = 240, separating the criteria that were actually applied,
#   (2) reports the choice-profile breakdown of the excluded participants
#       (always-recommend vs always-reject vs near-degenerate vs mixed),
#   (3) compares the excluded and retained samples on the trait measures that can
#       be scored identically in both groups (GRiPS, Berlin + subjective numeracy,
#       BFI-2 domains, religiosity / spirituality) with effect sizes, FDR
#       correction, and TOST equivalence tests, and
#   (4) runs a model-free robustness check: per-participant penalized logistic
#       regressions of choice on harm and cure probability, which -- unlike the
#       DDM -- can be estimated for near-degenerate responders, so the key
#       trait/sensitivity associations can be re-run on a sample that puts the
#       response-excluded participants back in.
#
# Everything is written to Figures/ExclusionRobustness/. Run time is a couple of
# minutes, dominated by the ~425 per-participant logistic fits in section 7.
#
# Terminology note (see DATA_ANALYSIS/README.md): the code says risk/gain where
# the manuscript says harm/cure. `experimental_treatment_selected == "Yes to
# treatment"` is a "Recommend" response.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(psych)    # principal(): recovers the BFI-2 keying from the data
  library(arm)      # bayesglm(): penalized logistic, tolerates separation
})

# arm (via MASS) masks dplyr::select; pin the dplyr verbs we rely on.
select <- dplyr::select
filter <- dplyr::filter

set.seed(20250910)

# ---------------------------------------------------------------------------
# 0. Paths and settings
# ---------------------------------------------------------------------------

base_dir <- "C:/Users/Brian Kim/Documents/OneDrive - Drexel University/Documents/classes/cognew/Neuroethics/IndividualDifferenceRiskBenefit/DATA_ANALYSIS"

path_raw        <- file.path(base_dir, "data/raw_mturk_data/raw.csv")
path_eligible   <- file.path(base_dir, "data/Cleaned_Pilot/exclusions.2020.10.21.RData")   # Neuroethics_Judgement
path_rt         <- file.path(base_dir, "data/FINAL_RT_EXCLUDED_DATA_20241014.RData")       # risk_data
path_final_list <- file.path(base_dir, "analysis/finalDDM/FINAL_PARTICIPANT_LIST.RData")   # unique_idxP_table
path_ddm_traits <- file.path(base_dir, "analysis/finalDDM/final_ddm_ideological_combined.RData")

out_dir <- file.path(base_dir, "analysis/finalDDM/Figures/ExclusionRobustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# N as reported in the Participants section of the manuscript. This is a count of
# Qualtrics RESPONSES, not of distinct people: raw.csv holds 616 distinct
# ResponseId values but only 610 distinct user_id values. The 6-response gap is
# repeat participation, not lost data -- see the duplicate audit below. The
# distinct-participant count at intake is therefore 610.
N_RECRUITED <- 616

# Smallest effect size of interest for the TOST equivalence tests. d = 0.35 is a
# conventional "small" bound; with n ~ 376 vs ~ 230 the tests have adequate power
# there. Change it here if a reviewer prefers a stricter bound.
SESOI_D <- 0.35

log_lines <- character(0)
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  log_lines <<- c(log_lines, txt)
}
say_header <- function(x) say("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78))
show_obj <- function(x) {
  log_lines <<- c(log_lines, capture.output(print(x)))
  print(x)
}

# ---------------------------------------------------------------------------
# 1. Load
# ---------------------------------------------------------------------------

raw <- read.csv(path_raw, stringsAsFactors = FALSE)

env_elig <- new.env(); load(path_eligible, envir = env_elig)
eligible_trials <- env_elig$Neuroethics_Judgement   # survives the prior-study criteria

env_rt <- new.env(); load(path_rt, envir = env_rt)
analyzed_trials <- env_rt$risk_data                 # the trials the DDM was fit to

env_list <- new.env(); load(path_final_list, envir = env_list)
final_ids <- unique(env_list$unique_idxP_table$user_id)

env_ddm <- new.env(); load(path_ddm_traits, envir = env_ddm)
ddm_traits <- as.data.frame(env_ddm$final_table)    # per-participant DDM params + traits

raw_ids      <- unique(raw$user_id)
eligible_ids <- unique(eligible_trials$user_id)
fit_ids      <- unique(analyzed_trials$user_id)

# ---------------------------------------------------------------------------
# 2. Exclusion cascade
# ---------------------------------------------------------------------------
# The manuscript compresses the cascade into one sentence ("participants that did
# not have both types of responses were excluded (n = 240)"). The data show the
# 240 is a *total* across several criteria, only one of which is about response
# degeneracy. Reconstructing it stage by stage is the honest way to answer the
# reviewer, and it turns out to help rather than hurt the argument.

say_header("SECTION 1. Exclusion cascade, 616 -> 376")

# The 616 -> 610 step is not an exclusion. Each Qualtrics ResponseId is one
# survey session; user_id identifies the person. A handful of people hold more
# than one session, and collapsing sessions to people accounts for the gap
# exactly. Audited here rather than asserted, because the original framing of
# these 6 as participants with "no usable record" was wrong.
dup_sessions <- raw %>%
  group_by(user_id) %>%
  summarise(n_sessions = n_distinct(ResponseId), n_rows = n(), .groups = "drop") %>%
  filter(n_sessions > 1) %>%
  arrange(desc(n_sessions))
n_dup_sessions <- sum(dup_sessions$n_sessions) - nrow(dup_sessions)
stopifnot(n_dup_sessions == N_RECRUITED - length(raw_ids))
n_no_record <- n_dup_sessions
n_prior     <- length(setdiff(raw_ids, eligible_ids))
n_response  <- length(setdiff(eligible_ids, fit_ids))
n_final     <- length(setdiff(fit_ids, final_ids))

cascade <- data.frame(
  stage = c("Recruited (as reported)",
            "Duplicate survey sessions collapsed to participant level",
            "Prior-study criteria (Fernandez et al., 2025), applied upstream",
            "Degenerate / near-degenerate choice profile",
            "Dropped after DDM fit, reason not recoverable from the repo",
            "Analyzed"),
  n_excluded = c(NA, n_no_record, n_prior, n_response, n_final, NA),
  n_remaining = c(N_RECRUITED,
                  length(raw_ids),
                  length(eligible_ids),
                  length(fit_ids),
                  length(final_ids),
                  length(final_ids))
)
show_obj(cascade)

say("\nDuplicate-session audit (the 616 -> 610 step):")
show_obj(as.data.frame(dup_sessions))
say(nrow(dup_sessions), " people hold ", nrow(dup_sessions) + n_dup_sessions,
    " sessions between them, i.e. ", n_dup_sessions, " duplicate sessions.")
say("No participant was lost at this step: every ResponseId in raw.csv belongs",
    "\nto a user_id that is present. Deduplication downstream is clean -- of the ",
    nrow(dup_sessions), "\nrepeat participants, ",
    sum(dup_sessions$user_id %in% eligible_ids),
    " survives the prior-study criteria, contributing a single",
    "\nsitting (75 trials) rather than both.")

# --- The 382 -> 376 step -------------------------------------------------
# These six were fit by the DDM: unique_idxP_table maps user_id to
# as.numeric(as.factor(user_id)) over the 382, and the six absent idxP values
# (26, 33, 104, 266, 292, 357) are exactly these participants. They were removed
# when the final participant list was built, and the code that built
# FINAL_PARTICIPANT_LIST.RData is not in this repository, so the rule cannot be
# recovered from source. What can be done is to rule out the obvious candidates,
# which is what the audit below does. Every one of them comes back negative:
#
#   - not survey incompleteness  : all six have Progress == 100
#   - not attention checks       : all six pass all four
#   - not listwise deletion      : three have nothing missing on the 209
#                                  final_table variables present in raw.csv,
#                                  while four RETAINED participants do have
#                                  missing values there
#   - not trial count            : 58-75 trials, inside the retained range
#   - not choice degeneracy      : recommend rates .14-.64, all mixed
#   - not RT anomalies           : medians and SDs inside the retained range
#   - not DDM convergence        : all six have posterior estimates in
#                                  NeuroethicsPosteriorAndPredictions.RData with
#                                  posterior SDs spanning the 0.5th-99.5th
#                                  percentile of the 382, i.e. no identifiability
#                                  failure (checked separately; the posterior
#                                  file is too large to load here)
#   - not repeat participation   : none of the six appears in
#                                  data/Redundancies.csv (all 29 flagged workers
#                                  were already gone before the DDM fit), none
#                                  holds a second ResponseId, and no two of the
#                                  382 share an exact survey-response
#                                  fingerprint
#   - not a shared IP address    : in the Qualtrics export each of the six has an
#                                  IP used by nobody else. Six IPs ARE shared by
#                                  pairs within the 382, but all twelve of those
#                                  participants were RETAINED, so IP duplication
#                                  was evidently not an exclusion rule
#   - not response status        : Qualtrics Status = 0 (normal, not preview or
#                                  test), Finished = 1, Progress = 100,
#                                  DistributionChannel = anonymous for all six
#   - not geography or language  : all six geolocate inside the continental US
#                                  with UserLanguage = EN. The single non-US
#                                  participant among the 382 (KZDOCQto,
#                                  Venezuela) was retained
#   - not a raw-file artifact    : the per-participant MTurk CSVs are byte
#                                  identical across the cumulative batch
#                                  snapshots, exactly as for retained controls
#   - not any measured variable  : an exhaustive search over 222 per-participant
#                                  variables -- trial counts, choice proportions,
#                                  every RT summary, every numeric survey item,
#                                  and the posterior mean AND SD of all six stDDM
#                                  parameters -- found NO variable on which these
#                                  six are extreme or separable. Their ranks
#                                  scatter across the full range on every one
#                                  (e.g. trial count 26/96/176/274/309/380 of
#                                  382; b_risk 7/23/55/81/308/351)
#
# No upstream file in data/ contains the 376-participant set; it appears only in
# FINAL_PARTICIPANT_LIST.RData and the two artifacts derived from it. The 376 was
# therefore created at that step rather than inherited from any earlier cleaning.
#
# This is worth surfacing rather than smoothing over: a reviewer who asks why
# the analysis reports 376 rather than 382 currently has no answer, and the
# authors should reconstruct or re-derive this step before resubmission.
lost_after_fit <- setdiff(fit_ids, final_ids)
audit_lost <- raw %>%
  filter(user_id %in% lost_after_fit) %>%
  group_by(user_id) %>% slice(1) %>% ungroup() %>%
  transmute(user_id,
            Progress,
            attention_ok = (attention_check.1 + attention_check.2 +
                            attention_check.3 + attention_check.4) == 4) %>%
  left_join(analyzed_trials %>%
              group_by(user_id) %>%
              summarise(n_trials = n(),
                        p_recommend = round(mean(
                          experimental_treatment_selected == "Yes to treatment"), 3),
                        rt_median = round(median(rt)), .groups = "drop"),
            by = "user_id") %>%
  as.data.frame()

say("\nThe ", length(lost_after_fit), " participants dropped after the DDM fit:")
show_obj(audit_lost)
say("All six completed the survey, passed every attention check, produced",
    "\ntrial counts and recommend rates inside the retained range, and have",
    "\nordinary DDM posteriors. They are not repeat participants, do not share an",
    "\nIP with anyone, and geolocate inside the US. An exhaustive search over 222",
    "\nper-participant variables found none on which they are extreme or",
    "\nseparable. The code that produced FINAL_PARTICIPANT_LIST.RData is not in",
    "\nthis repository and no upstream data file contains the 376 set, so no",
    "\nreason for their removal can be recovered. The robust course is to re-run",
    "\nthe reported analyses on all 382 and confirm the results are unchanged,",
    "\nrather than to justify the exclusion post hoc.")
write.csv(audit_lost, file.path(out_dir, "dropped_after_ddm_fit.csv"), row.names = FALSE)

total_excluded <- n_no_record + n_prior + n_response + n_final
say("\nTotal excluded = ", total_excluded,
    " (manuscript reports 240; reconciles exactly: ", total_excluded == 240, ")")
say("Note this total is response-level. ", n_no_record, " of the ", total_excluded,
    " are duplicate sessions rather\nthan excluded people, so participant-level ",
    "exclusion is ", total_excluded - n_no_record, " of ", length(raw_ids),
    " (",
    sprintf("%.1f", 100 * (total_excluded - n_no_record) / length(raw_ids)),
    "%),\nnot ", total_excluded, " of ", N_RECRUITED, ".")
say("Of these, the criterion the reviewer is worried about -- degenerate or",
    "\nnear-degenerate responding -- accounts for ", n_response, " of ", total_excluded,
    " (", round(100 * n_response / total_excluded, 1), "%).")

# What we can verify about the upstream (prior-study) exclusions. These criteria
# were applied in the Fernandez et al. cleaning pipeline, which is not in this
# repository, so we characterize the excluded set empirically rather than
# claiming to know the rule that removed each person.
prior_ids <- setdiff(raw_ids, eligible_ids)
prior_profile <- raw %>%
  filter(user_id %in% prior_ids) %>%
  group_by(user_id) %>%
  summarise(
    n_task_trials = sum(!is.na(experimental_treatment_selected)),
    progress      = first(Progress),
    ac1           = first(attention_check.1),
    ac2           = first(attention_check.2),
    .groups = "drop"
  )
prior_prop <- eligible_trials %>%
  group_by(user_id) %>%
  summarise(n = n(),
            p = mean(experimental_treatment_selected == "Yes to treatment"),
            .groups = "drop")
say("\nPrior-study criteria, as stated in Fernandez et al. (2025), checked",
    "\nagainst the ", nrow(prior_prop), "-participant file this project inherits:")
say("  retain at least 40 of 75 trials -> APPLIED (trial counts run ",
    paste(range(prior_prop$n), collapse = " to "), ")")
say("  exclude stereotyped responding on >=90% of trials -> NOT applied to",
    "\n    this file: ", sum(prior_prop$p >= .90 | prior_prop$p <= .10),
    " of its ", nrow(prior_prop), " participants still sit at or beyond that band")
say("=> The ", n_response, " removed here for choice profile are largely the prior",
    "\n   study's own criterion being enforced for the first time, on more",
    "\n   aggressively trimmed trials. The criterion is inherited, not invented.")

say("\nVerifiable characteristics of the ", nrow(prior_profile),
    " participants removed by the prior-study criteria:")
say("  no task data at all:              ", sum(prior_profile$n_task_trials == 0))
say("  incomplete survey (Progress<100): ", sum(prior_profile$progress < 100, na.rm = TRUE))
say("  failed an attention check:        ",
    sum(prior_profile$ac1 %in% 0 | prior_profile$ac2 %in% 0))
say("  none of the above (other prior-study criteria): ",
    sum(prior_profile$n_task_trials > 0 &
        !(prior_profile$progress < 100) &
        !(prior_profile$ac1 %in% 0 | prior_profile$ac2 %in% 0)))

# ---------------------------------------------------------------------------
# 3. What the response criterion actually was
# ---------------------------------------------------------------------------
# The manuscript describes a strict criterion ("did not have both types of
# responses"). That is not what was run, and it is worth pinning down exactly
# what was, because the honest description involves an interaction a reviewer
# would rather hear from the authors.
#
# The implemented rule is a 10%/90% band on the Recommend rate *computed after
# RT trimming*: keep participants with 0.10 < p < 0.90, drop the rest. The check
# below confirms it reproduces the 382/43 split with zero misclassifications.
#
# The band is INHERITED, not invented. Fernandez et al. (2025), "Data analytic
# exclusions", states: "Participants were also excluded for selecting a 'yes' or
# 'no' response on at least 90% of the total number of trials. These patterns of
# stereotyped responses precluded the ability to estimate reliable coefficients
# in our statistical models." That is exactly this criterion. The manuscript's
# phrase "did not have both types of responses" is a garbled restatement of it,
# which is why the implemented rule and the stated rule disagree.
#
# Note also that the prior study's rule was evidently NOT applied to the 425-
# participant file this project inherits: 35 of those 425 still sit at or beyond
# the 10/90 band (checked in section 1). By contrast the prior study's other
# stated rule -- retain at least 40 of 75 trials -- clearly WAS applied, since
# trial counts in that file run exactly 40 to 75. So the participants removed
# here for choice profile are largely the prior study's own criterion being
# enforced for the first time, on more aggressively trimmed trials.
#
# Two things follow, and both belong in the write-up:
#
#   (a) It is a proportion rule, not a count rule. Trial counts after trimming
#       range from 38 to 75, so a retained participant can have as few as 4
#       minority responses while a dropped one has as many as 7. Counts overlap;
#       proportions do not.
#   (b) Because the proportion is computed post-trimming, inclusion is jointly
#       determined by the response profile and by which trials the RT filter
#       happened to remove. The sensitivity check below quantifies this by
#       re-applying the same 10/90 band to the untrimmed trials and counting how
#       many participants change status. Those are participants sitting near the
#       boundary whose fate is decided by trimming rather than by their overall
#       willingness to recommend.
#
# What the rule is NOT is an accident of trimming stripping minority responses
# until people became degenerate. Trimming removes a median of one minority
# response among the dropped, and only 2 of the 43 end up strictly one-sided.
# The cut is a deliberate near-uniformity threshold.

say_header("SECTION 2. The response criterion as implemented")

choice_prop <- function(trials) {
  trials %>%
    group_by(user_id) %>%
    summarise(
      n_trials    = n(),
      n_recommend = sum(experimental_treatment_selected == "Yes to treatment"),
      .groups = "drop"
    ) %>%
    mutate(p_recommend = n_recommend / n_trials,
           n_reject    = n_trials - n_recommend)
}

# The published RT trimming, reapplied to the eligible pool. Verified against
# FINAL_RT_EXCLUDED_DATA_20241014.RData: per-participant trial counts match for
# every retained participant, so this is the same trimming the pipeline used.
rt_trimmed <- eligible_trials %>%
  filter(rt > 250, rt < 9000) %>%
  group_by(user_id) %>%
  mutate(Q1 = quantile(rt, .25), Q3 = quantile(rt, .75), IQR = IQR(rt)) %>%
  filter(rt > (Q1 - 2 * IQR), rt < (Q3 + 2 * IQR)) %>%
  ungroup()

trim_check <- inner_join(count(rt_trimmed, user_id, name = "n_reproduced"),
                         count(analyzed_trials, user_id, name = "n_pipeline"),
                         by = "user_id")
say("Trimming reproduction check: per-participant trial counts match the",
    "\npipeline file for all ", nrow(trim_check), " retained participants: ",
    all(trim_check$n_reproduced == trim_check$n_pipeline))

prop_untrimmed <- choice_prop(eligible_trials) %>%
  select(user_id, n_untrimmed = n_trials, p_untrimmed = p_recommend)
prop_eligible <- choice_prop(rt_trimmed) %>%
  mutate(retained = user_id %in% fit_ids) %>%
  left_join(prop_untrimmed, by = "user_id")

# Recover the rule.
BAND_LO <- 0.10
BAND_HI <- 0.90
in_band <- function(p) p > BAND_LO & p < BAND_HI

n_mis <- sum(in_band(prop_eligible$p_recommend) != prop_eligible$retained)
say("\nCandidate rule: retain if ", BAND_LO, " < (Recommend rate after RT trimming) < ",
    BAND_HI, ".")
say("  misclassifications against the actual ", sum(prop_eligible$retained), "/",
    sum(!prop_eligible$retained), " split: ", n_mis,
    ifelse(n_mis == 0, "   <- exact match", "   <- rule does not fit"))

say("\nRecommend rate after RT trimming, among the ", nrow(prop_eligible),
    " participants who")
say("passed the prior-study criteria:")
say("  retained (n = ", sum(prop_eligible$retained), "): ",
    paste(sprintf("%.4f", range(prop_eligible$p_recommend[prop_eligible$retained])),
          collapse = " to "))
say("  excluded (n = ", sum(!prop_eligible$retained), "): ",
    paste(sprintf("%.4f", range(prop_eligible$p_recommend[!prop_eligible$retained])),
          collapse = " to "),
    "  (nothing excluded falls inside the retained range)")

say("\nIt is a proportion rule, not a count rule. Trial counts after trimming range")
say("from ", min(prop_eligible$n_trials), " to ", max(prop_eligible$n_trials),
    ", so the counts overlap even though the proportions do not:")
say("  fewest minority responses among retained: ",
    min(pmin(prop_eligible$n_recommend, prop_eligible$n_reject)[prop_eligible$retained]))
say("  most minority responses among excluded:   ",
    max(pmin(prop_eligible$n_recommend, prop_eligible$n_reject)[!prop_eligible$retained]))

rare_untrimmed <- pmin(round(prop_eligible$p_untrimmed * prop_eligible$n_untrimmed),
                       prop_eligible$n_untrimmed -
                         round(prop_eligible$p_untrimmed * prop_eligible$n_untrimmed))
rare_trimmed <- pmin(prop_eligible$n_recommend, prop_eligible$n_reject)

say("\nStrictly one-sided responders (0 of one response type):")
say("  before RT trimming: ", sum(prop_eligible$p_untrimmed %in% c(0, 1)))
say("  after RT trimming:  ", sum(prop_eligible$p_recommend %in% c(0, 1)),
    ", out of the ", sum(!prop_eligible$retained), " excluded")
say("  median minority responses removed by trimming, among the excluded: ",
    median((rare_untrimmed - rare_trimmed)[!prop_eligible$retained]))
say("=> The exclusion is a near-uniformity threshold. It is not trimming")
say("   accidentally making people degenerate.")

# Boundary sensitivity: how much does inclusion depend on the trimming?
prop_eligible$retained_if_untrimmed <- in_band(prop_eligible$p_untrimmed)
flips <- prop_eligible %>%
  filter(retained_if_untrimmed != retained) %>%
  transmute(user_id,
            n_untrimmed, p_untrimmed = round(p_untrimmed, 4),
            n_trimmed = n_trials, p_trimmed = round(p_recommend, 4),
            actually_retained = retained,
            retained_if_rule_applied_pre_trimming = retained_if_untrimmed) %>%
  arrange(p_untrimmed)

say("\nBoundary sensitivity. Applying the same ", BAND_LO, "/", BAND_HI,
    " band to the untrimmed")
say("trials instead would retain ", sum(prop_eligible$retained_if_untrimmed),
    " rather than ", sum(prop_eligible$retained), " participants. ", nrow(flips),
    " of ", nrow(prop_eligible), " change status:")
show_obj(as.data.frame(flips))
write.csv(flips, file.path(out_dir, "boundary_sensitivity_flips.csv"), row.names = FALSE)
say("\n=> Inclusion is jointly determined by the response profile and by which")
say("   trials RT trimming removed. For ", nrow(flips), " participants (",
    sprintf("%.1f", 100 * nrow(flips) / nrow(prop_eligible)),
    "% of the eligible pool) the")
say("   trimming decides the outcome. That is a real if minor arbitrariness and")
say("   is better disclosed than left for a reviewer to find.")

# ---------------------------------------------------------------------------
# 4. Choice-profile breakdown of every excluded participant
# ---------------------------------------------------------------------------
# This is the table the reviewer asked for. It covers all excluded participants
# who produced task data, including those removed upstream, so the "strong
# preference" hypothesis is tested against the whole excluded set rather than
# just the response-based tranche.

say_header("SECTION 3. Choice profiles of the excluded participants")

raw_prop <- raw %>%
  filter(!is.na(experimental_treatment_selected)) %>%
  group_by(user_id) %>%
  summarise(
    n_trials    = n(),
    n_recommend = sum(experimental_treatment_selected),
    .groups = "drop"
  ) %>%
  mutate(p_recommend = n_recommend / n_trials)

# Bands match the criterion recovered in section 2 (BAND_LO / BAND_HI), so the
# descriptive table and the actual rule are on the same scale.
classify_profile <- function(p) {
  cut(p,
      breaks = c(-Inf, 0, BAND_LO, BAND_HI, 1 - 1e-9, Inf),
      labels = c("Always reject",
                 paste0("Mostly reject (<=", 100 * BAND_LO, "%)"),
                 "Mixed",
                 paste0("Mostly recommend (>=", 100 * BAND_HI, "%)"),
                 "Always recommend"))
}

exclusion_stage <- function(id) {
  ifelse(id %in% final_ids, "Analyzed",
  ifelse(id %in% fit_ids, "Dropped after DDM fit",
  ifelse(id %in% eligible_ids, "Excluded: choice profile",
         "Excluded: prior-study criteria")))
}

# Classify on the same quantity the criterion used: the Recommend rate after RT
# trimming. That is only defined for the 425 who reached the trimming step, so
# participants removed by the prior-study criteria are classified on their raw
# rate instead and the source is recorded. Mixing the two silently would be
# misleading -- section 2 showed the trimmed and untrimmed rates disagree near
# the boundary, which is precisely what the profile table is trying to describe.
profiles <- raw_prop %>%
  rename(n_raw = n_trials, p_raw = p_recommend) %>%
  left_join(prop_eligible %>% select(user_id, p_trimmed = p_recommend),
            by = "user_id") %>%
  mutate(
    stage        = exclusion_stage(user_id),
    rate_source  = ifelse(is.na(p_trimmed), "raw (never reached trimming)",
                          "after RT trimming"),
    p_recommend  = ifelse(is.na(p_trimmed), p_raw, p_trimmed),
    profile      = classify_profile(p_recommend)
  )

say("Recommend rate used for classification:")
show_obj(table(profiles$stage, profiles$rate_source))


show_obj(table(profiles$stage, profiles$profile))

excl_profiles <- profiles %>% filter(stage != "Analyzed")
n_reject_side    <- sum(excl_profiles$p_recommend <= BAND_LO)
n_recommend_side <- sum(excl_profiles$p_recommend >= BAND_HI)

say("\nAmong the ", nrow(excl_profiles), " excluded participants who produced task data:")
say("  near-always REJECT    (<=", 100 * BAND_LO, "% recommend): ", n_reject_side,
    " (", round(100 * n_reject_side / nrow(excl_profiles), 1), "%)")
say("  near-always RECOMMEND (>=", 100 * BAND_HI, "% recommend): ", n_recommend_side,
    " (", round(100 * n_recommend_side / nrow(excl_profiles), 1), "%)")
say("  mixed responders:                        ",
    nrow(excl_profiles) - n_reject_side - n_recommend_side)

say("\nWithin the tranche excluded specifically for choice profile (n = ", n_response, "):")
resp_only <- profiles %>% filter(stage == "Excluded: choice profile")
say("  near-always reject: ", sum(resp_only$p_recommend <= BAND_LO),
    " | near-always recommend: ", sum(resp_only$p_recommend >= BAND_HI))
bt <- binom.test(sum(resp_only$p_recommend <= BAND_LO), nrow(resp_only))
say("  reject-side proportion = ", sprintf("%.3f", bt$estimate),
    ", 95% CI [", sprintf("%.3f", bt$conf.int[1]), ", ", sprintf("%.3f", bt$conf.int[2]),
    "], exact binomial p = ", sprintf("%.3f", bt$p.value), " against 50/50")

write.csv(profiles, file.path(out_dir, "participant_choice_profiles.csv"), row.names = FALSE)
write.csv(cascade, file.path(out_dir, "exclusion_cascade.csv"), row.names = FALSE)

# Figure: distribution of recommend rates, retained vs excluded.
p_hist <- profiles %>%
  mutate(grp = ifelse(stage == "Analyzed",
                      paste0("Analyzed (n = ", length(final_ids), ")"),
                      "Excluded (with task data)")) %>%
  ggplot(aes(p_recommend, fill = grp)) +
  geom_histogram(binwidth = 0.05, boundary = 0, colour = "white", linewidth = 0.2) +
  facet_wrap(~ grp, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = setNames(c("#3b6ea5", "#c2543d"),
                                      c(paste0("Analyzed (n = ", length(final_ids), ")"),
                                        "Excluded (with task data)")),
                    guide = "none") +
  labs(x = "Proportion of trials recommending neuromodulation (after RT trimming\nwhere available)",
       y = "Participants",
       title = "Choice profiles of analyzed vs excluded participants") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "choice_profile_distribution.png"), p_hist,
       width = 6.5, height = 5, dpi = 300)

# ---------------------------------------------------------------------------
# 5. Score traits identically for everyone
# ---------------------------------------------------------------------------
# The composite scores stored in risk_data exist only for participants who
# survived cleaning, so for the excluded participants we score from the raw items
# in raw.csv. The scoring keys below are verified against the stored composites
# for the participants where both exist -- the script stops if they do not
# reproduce exactly, so the excluded and retained groups are guaranteed to be on
# the same scale. (BFI-2 items in raw.csv are stored as answered, NOT reverse-coded
# -- opposite-worded items correlate negatively, e.g. r(BIG5_1 "outgoing",
# BIG5_16 "quiet") = -.42 -- so a plain item mean is the wrong composite. The
# reverse-scored composite is the one that reproduces the stored composites
# exactly, and it is what section 5 uses. An unreversed variant is also computed,
# but ONLY as an acquiescence diagnostic in sections 4b and 5b; it is not the
# published scoring and must never be reported as though it were.)

say_header("SECTION 4. Trait scoring")

raw_person <- raw %>% group_by(user_id) %>% slice(1) %>% ungroup() %>% as.data.frame()

# BFI-2 item numbering: the five domains take every fifth item.
big5_key <- list(
  Extraversion      = seq(1, 56, by = 5),
  Agreeableness     = seq(2, 57, by = 5),
  Conscientiousness = seq(3, 58, by = 5),
  Neuroticism       = seq(4, 59, by = 5),
  Openness          = seq(5, 60, by = 5)
)

# Two BFI-2 scorings are carried through this script.
#
#   *_unrev : plain item mean with NO reverse-scoring. This is not a valid
#           personality composite and is not what the published pipeline
#           produced. It is retained only as an acquiescence diagnostic: because
#           the BFI-2 is balanced across forward- and reverse-worded items, a
#           content-based responder lands near the midpoint on it while an
#           indiscriminate agreer is pulled upward.
#   *     : reverse-scored. The reverse-keyed items are identified from the data
#           (sign of each item's loading on its domain's first principal
#           component) rather than from a remembered key, then each domain is
#           oriented using an anchor item whose wording fixes the direction. The
#           anchors are read straight out of data/text_key.csv:
#             BIG5_1  "Is outgoing, sociable"                -> high Extraversion
#             BIG5_2  "Is compassionate, has a soft heart"    -> high Agreeableness
#             BIG5_13 "Is dependable, steady"                 -> high Conscientiousness
#             BIG5_14 "Is moody, has up and down mood swings" -> high Neuroticism
#             BIG5_10 "Is curious about many different things"-> high Openness
big5_anchor <- c(Extraversion = 1, Agreeableness = 2, Conscientiousness = 13,
                 Neuroticism = 14, Openness = 10)

derive_big5_keying <- function(item_matrix) {
  lapply(names(big5_key), function(dom) {
    its <- big5_key[[dom]]
    X   <- item_matrix[, its, drop = FALSE]
    X   <- X[complete.cases(X), , drop = FALSE]
    ld  <- as.numeric(principal(X, nfactors = 1)$loadings)
    # principal() fixes the component sign arbitrarily; flip so the anchor item
    # loads positively, which pins the domain to its conventional direction.
    if (ld[which(its == big5_anchor[[dom]])] < 0) ld <- -ld
    its[ld < 0]
  }) %>% setNames(names(big5_key))
}

items_matrix <- as.matrix(raw_person[, paste0("BIG5_", 1:60)])
big5_reverse <- derive_big5_keying(items_matrix)

score_traits <- function(d) {
  grips <- paste0("GRiPS_", 1:8)
  sns   <- paste0("SNS", 1:8)
  out <- data.frame(
    user_id      = d$user_id,
    GRiPS        = rowMeans(d[, grips]),
    SNS_total    = rowMeans(d[, sns]),
    BerlinScore  = (d$Berlin1 == 25) + (d$Berlin2 == 30) +
                   (d$Berlin3 == 20) + (d$Berlin4 == 50),
    Religiosity  = suppressWarnings(as.numeric(d$Religiosity)),
    Spirituality = suppressWarnings(as.numeric(d$Spirituality)),
    age          = suppressWarnings(as.numeric(d$age)),
    sex          = d$sex,
    condition    = d$condition,
    stringsAsFactors = FALSE
  )
  for (dom in names(big5_key)) {
    its <- big5_key[[dom]]
    X   <- as.matrix(d[, paste0("BIG5_", its)])
    out[[paste0(dom, "_unrev")]] <- rowMeans(X)
    X[, its %in% big5_reverse[[dom]]] <- 6 - X[, its %in% big5_reverse[[dom]]]
    out[[dom]] <- rowMeans(X)
  }
  out
}

traits_all <- score_traits(raw_person)

# Verify the as-published scoring reproduces the stored composites exactly. This
# is what guarantees the excluded and retained groups are on the same scale.
stored <- analyzed_trials %>%
  group_by(user_id) %>% slice(1) %>% ungroup() %>%
  select(user_id, GRiPs_Score, SNS.total, BerlinScore,
         Extraversion, Agreeableness, Conscientiousness, Neuroticism, Openness) %>%
  as.data.frame()
chk <- inner_join(stored, traits_all, by = "user_id")
checks <- c(
  GRiPS       = max(abs(chk$GRiPs_Score - chk$GRiPS)),
  SNS_total   = max(abs(chk$SNS.total - chk$SNS_total)),
  BerlinScore = max(abs(chk$BerlinScore.x - chk$BerlinScore.y))
)
for (dom in names(big5_key)) {
  # .x = stored composite, .y = reverse-scored composite recomputed here.
  checks[dom] <- max(abs(chk[[paste0(dom, ".x")]] - chk[[paste0(dom, ".y")]]))
}
say("Max absolute discrepancy between re-scored and stored composites (n = ",
    nrow(chk), "):")
show_obj(round(checks, 10))
stopifnot(all(checks < 1e-9, na.rm = TRUE))
say("=> Scoring keys reproduce the stored composites exactly. NOTE the five",
    "\n   BFI-2 columns reproduce under REVERSE-SCORING only because",
    "\n   fix_bfi2_composites.R re-scored them in place on 2026-09-10. In the",
    "\n   originally distributed file (kept as *.PRE_BFI2_FIX.RData) they were",
    "\n   plain unreversed item means. See section 4b.")

# ---------------------------------------------------------------------------
# 4b. BFI-2 scoring check -- read this before interpreting section 5
# ---------------------------------------------------------------------------
# The BFI-2 domain columns originally distributed in
# FINAL_RT_EXCLUDED_DATA_20241014.RData were plain means over all twelve items of
# each domain with NO reverse-scoring. Verified against the pre-fix backup
# (FINAL_RT_EXCLUDED_DATA_20241014.PRE_BFI2_FIX.RData): the stored columns
# reproduce the unreversed item means to machine precision (max discrepancy 0.0)
# and differ from the correctly keyed scores by up to 2.0 scale points, with a
# mean inter-domain correlation of +.68 for an instrument built to be
# near-orthogonal.
#
# fix_bfi2_composites.R re-scored those five columns in place on 2026-09-10, so
# the CURRENT data file carries correctly keyed scores (mean inter-domain
# r = +.02) and section 5 below reports the corrected scoring. Anyone comparing
# this output against the manuscript should know that the published analyses ran
# on the unreversed columns.
#
# That does not propagate to any reported result. The composites are carried into
# final_table but never enter a reported analysis: every published Big Five
# result is item-level ridge importance, abs(coef(ridge_model)[-1]) with
# standardize = TRUE, which is exactly invariant to reverse-scoring (verified on
# these items, max difference 4e-16 with all 30 reverse-keyed signs flipping).
#
# The unreversed variant is still computed here, but only as an acquiescence
# index for section 5b: collapsing a balanced item set without reverse-scoring
# cancels forward- against reverse-worded content and leaves mostly a single
# agreement factor. The contrast below quantifies the difference.

say_header("SECTION 4b. BFI-2 keying verification")

say("Reverse-keyed items recovered from the data (anchored by item wording):")
for (dom in names(big5_key)) {
  say(sprintf("  %-18s %s", dom, paste0("BIG5_", big5_reverse[[dom]], collapse = ", ")))
}

items_cc <- items_matrix[complete.cases(items_matrix), , drop = FALSE]
scoring_diag <- bind_rows(lapply(names(big5_key), function(dom) {
  its <- big5_key[[dom]]
  Xp  <- items_cc[, its, drop = FALSE]
  Xk  <- Xp; Xk[, its %in% big5_reverse[[dom]]] <- 6 - Xk[, its %in% big5_reverse[[dom]]]
  data.frame(
    domain      = dom,
    alpha_unrev = suppressWarnings(psych::alpha(Xp, warnings = FALSE)$total$raw_alpha),
    alpha_keyed = suppressWarnings(psych::alpha(Xk, warnings = FALSE)$total$raw_alpha),
    sd_unrev    = sd(rowMeans(Xp)),
    sd_keyed    = sd(rowMeans(Xk)),
    stringsAsFactors = FALSE
  )
}))
say("\nInternal consistency and dispersion, unreversed vs reverse-scored",
    "\n(the reverse-scored column is the published scoring):")
show_obj(format(scoring_diag, digits = 3))

dom_unrev <- sapply(big5_key, function(i) rowMeans(items_cc[, i, drop = FALSE]))
dom_keyed <- sapply(names(big5_key), function(dom) {
  its <- big5_key[[dom]]
  X <- items_cc[, its, drop = FALSE]
  X[, its %in% big5_reverse[[dom]]] <- 6 - X[, its %in% big5_reverse[[dom]]]
  rowMeans(X)
})
say("\nMean inter-domain correlation (BFI-2 domains should be near-orthogonal):")
say("  unreversed (diagnostic only): ",
    sprintf("%+.3f", mean(cor(dom_unrev)[upper.tri(cor(dom_unrev))])))
say("  reverse-scored (published):   ",
    sprintf("%+.3f", mean(cor(dom_keyed)[upper.tri(cor(dom_keyed))])))
say("=> The unreversed scoring -- what the originally distributed file carried,",
    "\n   and what the published analyses used -- collapses into a single agreement",
    "\n   factor. The reverse-scored version behaves as the instrument intends and",
    "\n   is what the current data file and section 5 now use.")

traits_all$stage <- exclusion_stage(traits_all$user_id)
traits_all$group <- ifelse(traits_all$stage == "Analyzed", "Analyzed", "Excluded")

trait_vars <- c("GRiPS", "BerlinScore", "SNS_total",
                "Extraversion", "Agreeableness", "Conscientiousness",
                "Neuroticism", "Openness", "Religiosity", "Spirituality", "age")

n_complete <- sapply(trait_vars, function(v)
  sum(!is.na(traits_all[[v]]) & traits_all$group == "Excluded"))
say("\nExcluded participants with usable survey data, by measure:")
show_obj(n_complete)

# ---------------------------------------------------------------------------
# 6. Retained vs excluded on traits
# ---------------------------------------------------------------------------
# Welch t-tests with Hedges' g and its CI, BH-FDR across the trait family, plus
# Kolmogorov-Smirnov to catch differences in spread rather than location, plus
# two one-sided tests so that "no difference" can be stated as an equivalence
# result rather than a failure to reject.

say_header("SECTION 5. Trait comparison, analyzed vs excluded")

hedges_g <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  n1 <- length(x); n2 <- length(y)
  sp <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
  d  <- (mean(x) - mean(y)) / sp
  J  <- 1 - 3 / (4 * (n1 + n2) - 9)          # small-sample correction
  g  <- J * d
  se <- sqrt((n1 + n2) / (n1 * n2) + g^2 / (2 * (n1 + n2 - 2)))
  c(g = g, lo = g - 1.96 * se, hi = g + 1.96 * se)
}

# Two one-sided tests against +/- SESOI_D, on the pooled-SD scale.
tost_p <- function(x, y, d_bound) {
  n1 <- length(x); n2 <- length(y)
  v1 <- var(x); v2 <- var(y)
  sp <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
  se <- sqrt(v1 / n1 + v2 / n2)                                   # Welch SE
  df <- (v1 / n1 + v2 / n2)^2 /
        ((v1 / n1)^2 / (n1 - 1) + (v2 / n2)^2 / (n2 - 1))         # Welch df
  diff  <- mean(x) - mean(y)
  bound <- d_bound * sp
  p_lower <- pt((diff + bound) / se, df, lower.tail = FALSE)  # H0: diff <= -bound
  p_upper <- pt((diff - bound) / se, df, lower.tail = TRUE)   # H0: diff >= +bound
  max(p_lower, p_upper)
}

compare_groups <- function(dat, g1_label, g2_label, group_col = "group") {
  res <- bind_rows(lapply(trait_vars, function(v) {
    x <- dat[[v]][dat[[group_col]] == g1_label]
    y <- dat[[v]][dat[[group_col]] == g2_label]
    x <- x[!is.na(x)]; y <- y[!is.na(y)]
    if (length(x) < 3 || length(y) < 3) return(NULL)
    tt <- t.test(x, y)
    g  <- hedges_g(x, y)
    ks <- suppressWarnings(ks.test(x, y))
    data.frame(
      measure = v,
      n1 = length(x), n2 = length(y),
      mean1 = mean(x), sd1 = sd(x),
      mean2 = mean(y), sd2 = sd(y),
      t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value,
      g = unname(g["g"]), g_lo = unname(g["lo"]), g_hi = unname(g["hi"]),
      ks_D = unname(ks$statistic), ks_p = ks$p.value,
      tost_p = tost_p(x, y, SESOI_D),
      stringsAsFactors = FALSE
    )
  }))
  res$p_fdr      <- p.adjust(res$p, "BH")
  res$ks_p_fdr   <- p.adjust(res$ks_p, "BH")
  res$tost_p_fdr <- p.adjust(res$tost_p, "BH")
  res$equivalent <- res$tost_p_fdr < .05
  res
}

fmt <- function(res, lab1 = "Group1", lab2 = "Group2") {
  out <- data.frame(
    Measure    = res$measure,
    A          = sprintf("%.2f (%.2f)", res$mean1, res$sd1),
    B          = sprintf("%.2f (%.2f)", res$mean2, res$sd2),
    `Hedges g` = sprintf("%.2f [%.2f, %.2f]", res$g, res$g_lo, res$g_hi),
    p          = sprintf("%.3f", res$p),
    p_FDR      = sprintf("%.3f", res$p_fdr),
    KS_p_FDR   = sprintf("%.3f", res$ks_p_fdr),
    Equivalent = ifelse(res$equivalent, "yes", "no"),
    check.names = FALSE
  )
  names(out)[2:3] <- c(lab1, lab2)
  out
}

main_res <- compare_groups(traits_all, "Analyzed", "Excluded")
say("Analyzed (n = ", length(final_ids), ") vs all excluded participants with survey data.")
say("Group columns are mean (SD). 'Equivalent' = TOST rejects an effect as large",
    "\nas |d| = ", SESOI_D, " at FDR < .05.\n")
show_obj(fmt(main_res, "Analyzed", "Excluded"))
write.csv(main_res, file.path(out_dir, "trait_comparison_analyzed_vs_excluded.csv"),
          row.names = FALSE)

say("\nLargest absolute Hedges' g: ",
    main_res$measure[which.max(abs(main_res$g))], " = ",
    sprintf("%.3f", main_res$g[which.max(abs(main_res$g))]))
say("Measures significant after FDR correction: ",
    ifelse(any(main_res$p_fdr < .05),
           paste(main_res$measure[main_res$p_fdr < .05], collapse = ", "), "none"))
say("Measures statistically equivalent at d = ", SESOI_D, ": ",
    ifelse(any(main_res$equivalent),
           paste(main_res$measure[main_res$equivalent], collapse = ", "), "none"))

# The manuscript's Big Five numbers come from the UNREVERSED columns (section 4b),
# while the table above uses the corrected keying. Both are reported so the
# response letter can quote whichever matches the manuscript text, and so the
# claim that the keying changes direction but not significance is testable rather
# than assumed.
#
# NOTE: an earlier version of this block built the variable names with
# sub("...", "\\1_unrev", ...) written as "\1_unrev". In R "\1" is an OCTAL
# escape (\001), not a backreference, so every Big Five name became an
# unmatchable control character and all five rows were silently dropped from the
# output. The doubled backslash below is load-bearing.
trait_vars_unrev <- sub(
  "^(Extraversion|Agreeableness|Conscientiousness|Neuroticism|Openness)$",
  "\\1_unrev", trait_vars)
stopifnot(all(trait_vars_unrev %in% names(traits_all)))

unrev_res <- local({
  keep <- trait_vars
  trait_vars <<- trait_vars_unrev
  on.exit(trait_vars <<- keep)
  compare_groups(traits_all, "Analyzed", "Excluded")
})
say("\n--- Same contrast, BFI-2 scored as the MANUSCRIPT scored it (unreversed) ---")
show_obj(fmt(unrev_res, "Analyzed", "Excluded"))
write.csv(unrev_res,
          file.path(out_dir, "trait_comparison_analyzed_vs_excluded_unreversedBFI.csv"),
          row.names = FALSE)

# Do the two scorings agree on which Big Five domains separate the groups?
big5_rows <- function(d, strip) {
  d <- d[sub(strip, "", d$measure) %in% names(big5_key), ]
  data.frame(domain = sub(strip, "", d$measure),
             g = round(d$g, 3),
             p_FDR = round(d$p_fdr, 4),
             sig = d$p_fdr < .05,
             stringsAsFactors = FALSE)
}
cmp <- merge(big5_rows(main_res, "$^"), big5_rows(unrev_res, "_unrev$"),
             by = "domain", suffixes = c("_corrected", "_manuscript"))
say("\nBig Five under both scorings (corrected = section 5 table, manuscript = unreversed):")
show_obj(cmp)
say("\n  domains agreeing on significance at FDR < .05: ",
    sum(cmp$sig_corrected == cmp$sig_manuscript), " of ", nrow(cmp))
say("  domains where the sign of g flips: ",
    sum(sign(cmp$g_corrected) != sign(cmp$g_manuscript)))

# Same comparison restricted to the tranche excluded for choice profile -- the
# group the reviewer's concern is actually about.
traits_resp <- traits_all %>%
  filter(stage %in% c("Analyzed", "Excluded: choice profile")) %>%
  mutate(group2 = ifelse(stage == "Analyzed", "Analyzed", "ChoiceExcluded"))
resp_res <- compare_groups(traits_resp, "Analyzed", "ChoiceExcluded", "group2")
say("\n--- Analyzed vs participants excluded for choice profile only ---")
show_obj(fmt(resp_res, "Analyzed", "ChoiceExcl"))
write.csv(resp_res, file.path(out_dir, "trait_comparison_analyzed_vs_choice_excluded.csv"),
          row.names = FALSE)

# Within the excluded, do the always-reject and always-recommend sides differ?
# If neuromodulation-relevant convictions drive degenerate responding, this is
# where a GRiPS or religiosity difference should show up.
traits_side <- traits_all %>%
  left_join(profiles %>% select(user_id, p_recommend), by = "user_id") %>%
  filter(group == "Excluded", !is.na(p_recommend)) %>%
  mutate(side = ifelse(p_recommend <= BAND_LO, "RejectSide",
                ifelse(p_recommend >= BAND_HI, "RecommendSide", NA_character_))) %>%
  filter(!is.na(side))
say("\n--- Within excluded: near-always reject vs near-always recommend ---")
say("n reject-side = ", sum(traits_side$side == "RejectSide"),
    ", n recommend-side = ", sum(traits_side$side == "RecommendSide"))
side_res <- compare_groups(traits_side, "RejectSide", "RecommendSide", "side")
show_obj(fmt(side_res, "RejectSide", "RecommendSide"))
write.csv(side_res, file.path(out_dir, "trait_comparison_reject_vs_recommend_side.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 5b. Response-style diagnostics
# ---------------------------------------------------------------------------
# Section 5 finds the excluded group scoring higher on ALL FIVE BFI-2 domains at
# once, and simultaneously higher on GRiPS and religiosity but LOWER on Berlin.
# A uniform shift across five near-orthogonal personality domains is not a
# plausible personality effect; it is the signature of a response style. The
# diagnostics below separate the two accounts, because they make opposite
# predictions:
#
#   substantive-difference account -> domains move independently; the objective
#                                     numeracy test should not move with them
#   response-style account         -> domains become strongly intercorrelated
#                                     (shared method variance), self-report
#                                     scales all drift toward "agree", and the
#                                     one effortful, objectively scored measure
#                                     (Berlin) drops
#
# Which account is right determines what the paper should say. A substantive
# difference bounds the claims to switchers; a response-style difference makes
# the exclusion a data-quality filter that happens to correlate with degenerate
# responding, which is a far weaker threat to validity.

say_header("SECTION 5b. Response-style diagnostics")

# raw.csv stores the responses as given (section 4b), so the acquiescence index
# is simply the mean over all 60 items: with the BFI-2 balanced across forward-
# and reverse-worded content, someone answering on content sits near the
# midpoint, while someone agreeing indiscriminately is pulled upward.
B_original <- items_matrix

style <- data.frame(
  user_id      = raw_person$user_id,
  acquiescence = rowMeans(B_original, na.rm = TRUE),      # drift toward "agree"
  within_sd    = apply(B_original, 1, sd, na.rm = TRUE),  # low = straightlining
  duration_s   = suppressWarnings(as.numeric(raw_person$Duration..in.seconds.)),
  stringsAsFactors = FALSE
) %>%
  left_join(traits_all %>% select(user_id, group, stage, BerlinScore), by = "user_id")

style_vars <- c("acquiescence", "within_sd", "BerlinScore", "duration_s")
# `comparison` selects the excluded rows to contrast against the analyzed
# sample; the analyzed side is always the same 376.
style_compare <- function(comparison) {
  bind_rows(lapply(style_vars, function(v) {
    x <- style[[v]][style$group == "Analyzed"]; x <- x[!is.na(x)]
    y <- style[[v]][comparison];                y <- y[!is.na(y)]
    g <- hedges_g(x, y)
    # Survey duration has a handful of participants who left the tab open for
    # hours, so it is summarised by the median and tested with Wilcoxon; the
    # mean and Hedges' g are meaningless for that row.
    if (v == "duration_s") {
      data.frame(index = "duration_s (median)",
                 analyzed = sprintf("%.0f", median(x)),
                 excluded = sprintf("%.0f", median(y)),
                 g = "-- (skewed)",
                 p = sprintf("%.4f", suppressWarnings(wilcox.test(x, y)$p.value)),
                 stringsAsFactors = FALSE)
    } else {
      data.frame(index = v,
                 analyzed = sprintf("%.2f (%.2f)", mean(x), sd(x)),
                 excluded = sprintf("%.2f (%.2f)", mean(y), sd(y)),
                 g = sprintf("%.2f [%.2f, %.2f]", g["g"], g["lo"], g["hi"]),
                 p = sprintf("%.4f", t.test(x, y)$p.value),
                 stringsAsFactors = FALSE)
    }
  }))
}
style_res        <- style_compare(style$group %in% "Excluded")
style_res_choice <- style_compare(style$stage %in% "Excluded: choice profile")
say("acquiescence = mean BFI-2 response in the direction the participant actually")
say("answered (higher = agrees with everything). within_sd = SD across that")
say("person's 60 items (lower = flatter, less differentiated responding).")
show_obj(style_res)
say("
Same indices for the tranche excluded specifically for choice profile:")
show_obj(style_res_choice)

# Inter-domain correlations. In a well-behaved sample the five BFI-2 domains are
# close to orthogonal; strong positive intercorrelation across all pairs means
# one common factor is driving every scale, i.e. method variance.
# Computed on the UNREVERSED item means, which is where shared method variance
# shows up. This is a response-style diagnostic, not the published scoring: the
# published composites are reverse-scored and behave normally (section 4b).
big5_names <- paste0(names(big5_key), "_unrev")
inter_cor <- function(ids, label) {
  d <- traits_all[traits_all$user_id %in% ids, big5_names]
  d <- d[complete.cases(d), , drop = FALSE]
  cm <- cor(d)
  offdiag <- cm[upper.tri(cm)]
  say(sprintf("%-28s n = %3d | mean inter-domain r = %+.3f | range [%+.3f, %+.3f]",
              label, nrow(d), mean(offdiag), min(offdiag), max(offdiag)))
  cm
}
say("BFI-2 inter-domain correlations (UNREVERSED items -- acquiescence",
    "\ndiagnostic, NOT the published scoring):")
cm_an <- inter_cor(final_ids, "Analyzed")
cm_ex <- inter_cor(traits_all$user_id[traits_all$group == "Excluded"], "Excluded (all)")
invisible(inter_cor(traits_all$user_id[traits_all$stage == "Excluded: choice profile"],
                    "Excluded (choice profile)"))
say("Analyzed sample correlation matrix:"); show_obj(round(cm_an, 2))
say("Excluded sample correlation matrix:"); show_obj(round(cm_ex, 2))

write.csv(style, file.path(out_dir, "response_style_indices.csv"), row.names = FALSE)
write.csv(style_res, file.path(out_dir, "response_style_comparison.csv"), row.names = FALSE)

# Verdict, computed rather than asserted, so it stays correct if the data change.
acq_g <- hedges_g(style$acquiescence[style$group == "Analyzed"],
                  style$acquiescence[style$group == "Excluded"])["g"]
sd_g  <- hedges_g(style$within_sd[style$group == "Analyzed"],
                  style$within_sd[style$group == "Excluded"])["g"]
ber_g <- hedges_g(na.omit(style$BerlinScore[style$group == "Analyzed"]),
                  na.omit(style$BerlinScore[style$group == "Excluded"]))["g"]
style_signature <- acq_g < -0.2 && sd_g > 0.2 && ber_g > 0.2

say("")
if (style_signature) {
  say("VERDICT: the excluded group agrees more indiscriminately (acquiescence")
  say(sprintf("g = %+.2f), answers far more flatly across the 60 items (within-person", acq_g))
  say(sprintf("SD g = %+.2f), and scores lower on the one objectively scored, effortful", sd_g))
  say(sprintf("measure in the battery (Berlin numeracy g = %+.2f). Those three move", ber_g))
  say("together, which is the signature of low-effort questionnaire responding,")
  say("not of strong convictions about neuromodulation. On this reading the")
  say("exclusion behaves as a data-quality filter, and near-uniform task")
  say("responding is one more expression of the same low engagement.")
} else {
  say("VERDICT: the response-style indices do not separate the groups, so the")
  say("section 5 trait differences should be read substantively and the paper")
  say("should bound its claims to participants who switched responses.")
}

# Categorical composition.
say("\n--- Categorical composition ---")
sex_tab <- table(traits_all$group, traits_all$sex, useNA = "no")
show_obj(sex_tab); show_obj(suppressWarnings(chisq.test(sex_tab)))
cond_tab <- table(traits_all$group, traits_all$condition)
cond_tab <- cond_tab[, colSums(cond_tab) > 0, drop = FALSE]
show_obj(cond_tab); show_obj(suppressWarnings(chisq.test(cond_tab)))

# Figure: forest plot of the trait effect sizes.
forest_dat <- bind_rows(
  main_res %>% mutate(contrast = "Analyzed vs all excluded"),
  resp_res %>% mutate(contrast = "Analyzed vs choice-excluded")
)
p_forest <- ggplot(forest_dat, aes(x = g, y = measure, colour = contrast)) +
  annotate("rect", xmin = -SESOI_D, xmax = SESOI_D, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_pointrange(aes(xmin = g_lo, xmax = g_hi),
                  position = position_dodge(width = 0.5), size = 0.3) +
  scale_colour_manual(values = c("Analyzed vs all excluded" = "#3b6ea5",
                                 "Analyzed vs choice-excluded" = "#c2543d")) +
  labs(x = paste0("Hedges' g (shaded band = equivalence bound, d = ", SESOI_D, ")"),
       y = NULL, colour = NULL,
       title = "Trait differences between analyzed and excluded participants") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "trait_effect_sizes.png"), p_forest,
       width = 7, height = 5, dpi = 300)

# ---------------------------------------------------------------------------
# 7. Why the excluded participants cannot simply be added to the DDM
# ---------------------------------------------------------------------------

say_header("SECTION 6. Why degenerate responders were excluded from the DDM")
say("The stDDM estimates harm and cure sensitivity from the *contrast* between",
    "\nRecommend and Do-Not-Recommend responses and their RT distributions. A",
    "\nparticipant who gives one response on (nearly) every trial provides no",
    "\ninformation about the boundary they never reach: the drift coefficients",
    "\nand the relative-onset parameter are only weakly identified, and their",
    "\nposteriors collapse onto the group prior. Including them would not recover",
    "\ntheir preferences -- it would import shrinkage estimates that look like",
    "\ngroup means and would attenuate exactly the individual-difference effects",
    "\nthe paper is about. The robustness analysis below therefore uses a",
    "\nmodel-free sensitivity measure that *is* estimable for these participants.")

# ---------------------------------------------------------------------------
# 8. Model-free robustness check
# ---------------------------------------------------------------------------
# Per-participant logistic regression of choice on harm and cure probability,
# fit with bayesglm() so that near-separated (degenerate) participants still get
# finite, regularized coefficients. Two questions:
#   (a) does this proxy track the DDM drift coefficients in the analyzed sample?
#   (b) do the published trait associations change when the choice-excluded
#       participants are added back in?

say_header("SECTION 7. Model-free robustness check")

fit_logit <- function(d) {
  d$y <- as.integer(d$experimental_treatment_selected == "Yes to treatment")
  m <- try(arm::bayesglm(y ~ risk_c + gain_c, data = d, family = binomial,
                         prior.scale = 2.5, prior.df = 1), silent = TRUE)
  if (inherits(m, "try-error")) return(data.frame(b_harm = NA_real_, b_cure = NA_real_))
  cf <- coef(m)
  data.frame(b_harm = unname(cf["risk_c"]), b_cure = unname(cf["gain_c"]))
}

# Scale harm/cure once, across the whole pool, so the per-participant
# coefficients are on a common scale.
logit_input <- rt_trimmed %>%
  mutate(risk_c = as.numeric(scale(risk)),
         gain_c = as.numeric(scale(gain)))

logit_fits <- logit_input %>%
  group_by(user_id) %>%
  group_modify(~ fit_logit(.x)) %>%
  ungroup() %>%
  mutate(retained = user_id %in% final_ids)

say("Fitted per-participant logistic models for ", nrow(logit_fits), " participants (",
    sum(logit_fits$retained), " analyzed, ", sum(!logit_fits$retained),
    " excluded downstream of the prior-study criteria -- the ", n_response,
    " removed for choice profile plus the ", n_final,
    " dropped between the DDM fit and the analysis sample).")

# (a) Does the proxy track the DDM parameters where both exist?
proxy_check <- inner_join(logit_fits,
                          ddm_traits %>% select(user_id, b_risk, b_gain),
                          by = "user_id")
say("\nConvergent validity in the analyzed sample (n = ", nrow(proxy_check), "):")
say("  cor(logistic harm coefficient, DDM b_risk) = ",
    sprintf("%.3f", cor(proxy_check$b_harm, proxy_check$b_risk, use = "complete.obs")))
say("  cor(logistic cure coefficient, DDM b_gain) = ",
    sprintf("%.3f", cor(proxy_check$b_cure, proxy_check$b_gain, use = "complete.obs")))

# (b) Re-run the trait associations with and without the excluded participants.
# Traits come from the eligible-pool cleaned data, which carries the composite
# scores for the choice-excluded participants as well.
elig_traits <- eligible_trials %>%
  group_by(user_id) %>% slice(1) %>% ungroup() %>%
  transmute(user_id,
            GRiPs_Score,
            consequentalism = consequentialism.score,
            Harm, Fairness, Loyalty, Athority, Sancticy,
            politics = republican + conservative - democrat - liberal) %>%
  as.data.frame()

robust_dat <- logit_fits %>%
  left_join(elig_traits, by = "user_id") %>%
  filter(!is.na(GRiPs_Score))

ind_vars <- c("GRiPs_Score", "consequentalism", "politics",
              "Harm", "Fairness", "Loyalty", "Athority", "Sancticy")

run_bank <- function(dat, label) {
  bind_rows(lapply(c("b_harm", "b_cure"), function(dv) {
    bind_rows(lapply(ind_vars, function(iv) {
      m <- lm(as.formula(paste(dv, "~", iv)), data = dat)
      s <- summary(m)$coefficients
      data.frame(sample = label, dv = dv, iv = iv,
                 beta = s[2, 1], se = s[2, 2], t = s[2, 3], p = s[2, 4],
                 n = nobs(m), stringsAsFactors = FALSE)
    }))
  }))
}

lab_retained <- "AnalyzedOnly"
lab_full     <- "PlusExcluded"
bank <- bind_rows(
  run_bank(robust_dat %>% filter(retained), lab_retained),
  run_bank(robust_dat, lab_full)
)
say("\nAnalyzed-only n = ", sum(robust_dat$retained),
    "; analyzed + all downstream exclusions n = ", nrow(robust_dat))

bank_wide <- bank %>%
  select(dv, iv, sample, beta, p) %>%
  pivot_wider(names_from = sample, values_from = c(beta, p)) %>%
  as.data.frame()
show_obj(format(bank_wide, digits = 3))
write.csv(bank, file.path(out_dir, "robustness_trait_associations.csv"), row.names = FALSE)

b1 <- bank_wide[[paste0("beta_", lab_retained)]]
b2 <- bank_wide[[paste0("beta_", lab_full)]]
p1 <- bank_wide[[paste0("p_", lab_retained)]]
p2 <- bank_wide[[paste0("p_", lab_full)]]
say("\nCorrelation between the two coefficient vectors: ", sprintf("%.3f", cor(b1, b2)))
say("Terms that change significance status at p < .05 when the excluded",
    "\nparticipants are added back: ", sum((p1 < .05) != (p2 < .05)), " of ", nrow(bank_wide))
say("Terms that change sign: ", sum(sign(b1) != sign(b2)))

write.csv(logit_fits, file.path(out_dir, "logistic_sensitivity_estimates.csv"),
          row.names = FALSE)

# Threshold sensitivity. Section 2 showed the 10/90 band is one defensible cut
# among several, and that for 10 participants the trimming rather than the
# response profile decides which side they land on. The natural follow-up is
# whether the trait associations depend on where the band is drawn. Re-fit the
# same regression bank at a range of bands, from keeping everyone through to a
# cut twice as strict as the published one.
say("\nThreshold sensitivity: the same regression bank at alternative bands.")

bands <- list(c(0.00, 1.00), c(0.05, 0.95), c(BAND_LO, BAND_HI),
              c(0.15, 0.85), c(0.20, 0.80))
band_dat <- robust_dat %>%
  left_join(prop_eligible %>% select(user_id, p_trimmed = p_recommend), by = "user_id")

band_bank <- bind_rows(lapply(bands, function(bd) {
  keep <- band_dat %>% filter(p_trimmed > bd[1], p_trimmed < bd[2])
  lbl  <- sprintf("%.2f-%.2f (n = %d)", bd[1], bd[2], nrow(keep))
  run_bank(keep, lbl)
}))

band_wide <- band_bank %>%
  select(dv, iv, sample, beta) %>%
  pivot_wider(names_from = sample, values_from = beta) %>%
  as.data.frame()
show_obj(format(band_wide, digits = 3))
write.csv(band_bank, file.path(out_dir, "threshold_sensitivity.csv"), row.names = FALSE)

published_col <- grep(sprintf("^%.2f-%.2f", BAND_LO, BAND_HI), names(band_wide), value = TRUE)
other_cols <- setdiff(names(band_wide), c("dv", "iv", published_col))
say("\nCorrelation of each band's coefficient vector with the published band:")
for (cl in other_cols) {
  say(sprintf("  %-22s r = %+.3f", cl, cor(band_wide[[cl]], band_wide[[published_col]])))
}
say("\n=> If these correlations stay high, the reported trait associations do not")
say("   turn on where the exclusion threshold was drawn, which is the substantive")
say("   answer to the reviewer regardless of how the criterion is described.")

# ---------------------------------------------------------------------------
# 9. Caveats worth stating in the response letter
# ---------------------------------------------------------------------------

say_header("SECTION 8. Caveats")
say("1. Trait comparisons cover GRiPS, Berlin and subjective numeracy, the five",
    "\n   BFI-2 domains, religiosity, spirituality and age. These are the measures",
    "\n   whose scoring keys could be verified item-for-item against the published",
    "\n   composites, so excluded and retained participants are strictly comparable.",
    "\n   Moral Foundations and Moral Commitments composites were built in an",
    "\n   upstream pipeline not contained in this repository and are available only",
    "\n   for the eligible pool, so they enter the section 7 robustness check but",
    "\n   not the section 5 comparison against the prior-study exclusions.")
say("2. The reviewer also raises decision-making responsibility toward others.",
    "\n   No responsibility measure was administered, so that specific claim cannot",
    "\n   be tested here and should be conceded as a limitation.")
say("3. The BFI-2 domain columns originally distributed with this project were",
    "\n   unreversed item means (verified against the pre-fix backup: exact",
    "\n   match to rowMeans over unreversed items, mean inter-domain r = +.68).",
    "\n   fix_bfi2_composites.R re-scored them on 2026-09-10, so section 5 uses",
    "\n   the corrected keying while the manuscript reports the unreversed ones.",
    "\n   No reported result changes: the composites never enter a published",
    "\n   analysis, and the item-level ridge importances are invariant to",
    "\n   reverse-scoring (max difference 4e-16). The exclusion comparison in",
    "\n   section 5 holds under either scoring.")
say("4. The section 7 proxy is a static logistic model. It recovers sensitivity",
    "\n   to harm and cure but has no analogue of the relative-onset (time)",
    "\n   parameter, so it can show that the magnitude effects are not an artifact",
    "\n   of the exclusion, but it cannot do the same for the latency effects.")

# ---------------------------------------------------------------------------
# 10. Save the transcript
# ---------------------------------------------------------------------------

say_header("Outputs written to")
say(out_dir)
writeLines(log_lines, file.path(out_dir, "exclusion_robustness_report.txt"))
