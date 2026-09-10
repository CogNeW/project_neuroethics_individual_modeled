# ---------------------------------------------------------------------------
# Re-score the BFI-2 domain composites in FINAL_RT_EXCLUDED_DATA_20241014.RData
# ---------------------------------------------------------------------------
# The five BFI-2 domain columns in `risk_data` (Extraversion, Agreeableness,
# Conscientiousness, Neuroticism, Openness) were stored as plain means over all
# twelve items of each domain, with no reverse-scoring. Because each BFI-2 domain
# is balanced across forward- and reverse-worded items, the reverse-worded items
# cancel the forward-worded ones and what remains is largely acquiescence: the
# five domains correlate at mean r = +.80 for an instrument built to be
# near-orthogonal, alpha falls to .62-.70, and the SDs are compressed ~24%.
#
# Evidence the stored columns are unreversed means:
#   - items in data/raw_mturk_data/raw.csv are raw 1-5 responses; opposite-worded
#     items correlate negatively, e.g. r(BIG5_1 "Is outgoing, sociable",
#     BIG5_16 "Tends to be quiet") = -.42
#   - the stored columns reproduce rowMeans() over the unreversed items to
#     machine precision for all 382 analyzed participants (asserted below)
#
# This script replaces those five columns with correctly keyed scores. It does
# NOT change any manuscript result: the composites are carried into final_table
# but never enter a reported analysis, and every published Big Five result is
# item-level ridge importance -- abs(coef()) with standardize = TRUE -- which is
# exactly invariant to reverse-scoring. The fix is so that anyone reusing the
# distributed data does not inherit mis-scored domain scores.
#
# The original file is copied to *.PRE_BFI2_FIX.RData before anything is written.
# Run from the repository root.

suppressMessages(library(dplyr))

data_path <- "data/FINAL_RT_EXCLUDED_DATA_20241014.RData"
raw_path  <- "data/raw_mturk_data/raw.csv"
stopifnot(file.exists(data_path), file.exists(raw_path))

# BFI-2 domains: every fifth item. Reverse-keyed items per Soto & John (2017);
# each domain has exactly 6 of its 12 items reverse-worded.
big5_key <- list(
  Extraversion      = seq(1, 56, by = 5),
  Agreeableness     = seq(2, 57, by = 5),
  Conscientiousness = seq(3, 58, by = 5),
  Neuroticism       = seq(4, 59, by = 5),
  Openness          = seq(5, 60, by = 5)
)
big5_reverse <- c(3, 4, 5, 8, 9, 11, 12, 16, 17, 22, 23, 24, 25, 26, 28, 29, 30,
                  31, 36, 37, 42, 44, 45, 47, 48, 49, 50, 51, 55, 58)
stopifnot(all(sapply(big5_key, function(it) sum(it %in% big5_reverse)) == 6))

raw <- read.csv(raw_path, stringsAsFactors = FALSE)
rp  <- raw[!duplicated(raw$user_id), ]
items <- sapply(rp[, paste0("BIG5_", 1:60)], function(z) suppressWarnings(as.numeric(z)))
rownames(items) <- rp$user_id
stopifnot(all(items >= 1 & items <= 5, na.rm = TRUE))   # 1-5 Likert, as expected

load(data_path)                                          # -> risk_data
stopifnot(exists("risk_data"))
idx <- match(risk_data$user_id, rownames(items))
stopifnot(!any(is.na(idx)))

cat("Verifying the stored columns are unreversed item means...\n")
for (dom in names(big5_key)) {
  unrev <- rowMeans(items[idx, big5_key[[dom]], drop = FALSE])
  d <- max(abs(risk_data[[dom]] - unrev), na.rm = TRUE)
  cat(sprintf("  %-18s max|stored - unreversed| = %.15f\n", dom, d))
  stopifnot(d < 1e-9)
}
cat("=> confirmed.\n\n")

backup <- sub("[.]RData$", ".PRE_BFI2_FIX.RData", data_path)
if (!file.exists(backup)) {
  file.copy(data_path, backup)
  cat("Original copied to", backup, "\n\n")
}

cat("Re-scoring with the published BFI-2 key...\n")
before <- sapply(names(big5_key), function(d) risk_data[[d]])
for (dom in names(big5_key)) {
  it <- big5_key[[dom]]
  X  <- items[idx, it, drop = FALSE]
  fl <- it %in% big5_reverse
  X[, fl] <- 6 - X[, fl]                                 # 1-5 scale: reverse is 6 - x
  risk_data[[dom]] <- rowMeans(X)
}
after <- sapply(names(big5_key), function(d) risk_data[[d]])

per_person <- !duplicated(risk_data$user_id)
cb <- cor(before[per_person, ]); ca <- cor(after[per_person, ])
cat(sprintf("\n  mean inter-domain r : %+.3f (before) -> %+.3f (after)\n",
            mean(cb[upper.tri(cb)]), mean(ca[upper.tri(ca)])))
cat(sprintf("  mean domain SD      :  %.3f (before) ->  %.3f (after)\n",
            mean(apply(before[per_person, ], 2, sd)),
            mean(apply(after[per_person, ], 2, sd))))

save(risk_data, file = data_path)
cat("\nWrote", data_path, "\n")
