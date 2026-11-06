#!/usr/bin/env Rscript
# Compare the subject-weight structure of the dec2024 and oct2026 samples, to explain why
# scen3's Step 2 (sub_s1) estimates var_subid without bias on dec2024 but not on oct2026.

repo_root <- "/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations"
data_dirs <- c(dec2024 = "/work/users/b/e/beibo/v3data/Dec2024/sample_miss00pct",
               oct2026 = "/work/users/b/e/beibo/v3data/Oct2026/sample_miss00pct")
sim_ids <- 1:50
keep_cols <- c("BGid", "hhid", "subid", "age_base", "y_gfr", "W_bghh", "W_sub", "sub_s1", "sub_s2")

sample_stats <- function(file) {
  visits <- read.csv(file)[keep_cols]
  visits$hh  <- paste(visits$BGid, visits$hhid)
  visits$sub <- paste(visits$hh, visits$subid)
  subj <- visits[!duplicated(visits$sub), c("hh", "sub", "age_base", "W_bghh", "W_sub", "sub_s1", "sub_s2")]
  subj$y_mean <- tapply(visits$y_gfr, visits$sub, mean)[subj$sub]

  n_j <- ave(subj$W_sub, subj$hh, FUN = length)
  multi <- n_j > 1
  hh_cv <- tapply(subj$W_sub[multi], subj$hh[multi], function(w) sd(w) / mean(w))
  # Within-household informativeness: subject weight vs. subject deviation from household mean outcome
  y_dev <- subj$y_mean - ave(subj$y_mean, subj$hh)
  lw_dev <- log(subj$W_sub) - ave(log(subj$W_sub), subj$hh)

  c(n_hh = length(unique(subj$hh)),
    n_sub = nrow(subj),
    sub_per_hh = nrow(subj) / length(unique(subj$hh)),
    pct_hh_single = mean(tapply(subj$sub, subj$hh, length) == 1),
    pct_sub_in_multi_hh = mean(multi),
    pct_multi_hh_unequal_wsub = mean(hh_cv > 1e-8),
    mean_cv_wsub_within_hh = mean(hh_cv),
    cv_wsub_overall = sd(subj$W_sub) / mean(subj$W_sub),
    cv_wbghh_overall = sd(subj$W_bghh) / mean(subj$W_bghh),
    mean_sub_s1_multi = mean(subj$sub_s1[multi]),
    mean_sub_s2_multi = mean(subj$sub_s2[multi]),
    mean_abs_s1_minus_s2 = mean(abs(subj$sub_s1 - subj$sub_s2)[multi]),
    wsub_ratio_age45_vs_young = mean(subj$W_sub[subj$age_base >= 45]) / mean(subj$W_sub[subj$age_base < 45]),
    cor_within_hh_logwsub_ydev = cor(lw_dev[multi], y_dev[multi]),
    cor_overall_logwsub_y = cor(log(subj$W_sub), subj$y_mean))
}

summary_rows <- lapply(names(data_dirs), function(pop) {
  stats <- sapply(sim_ids, function(i) sample_stats(file.path(data_dirs[[pop]], sprintf("samplemiss_%d.csv", i))))
  data.frame(population = pop, metric = rownames(stats), mean = rowMeans(stats),
             sd_across_samples = apply(stats, 1, sd), n_samples = ncol(stats))
})
out <- do.call(rbind, summary_rows)
wide <- reshape(out[c("population", "metric", "mean")], idvar = "metric", timevar = "population", direction = "wide")

dir.create(file.path(repo_root, "data", "summaries", "diagnostics"), showWarnings = FALSE)
write.csv(out, file.path(repo_root, "data", "summaries", "diagnostics", "weight_comparison_dec2024_oct2026.csv"),
          row.names = FALSE)
print(wide, row.names = FALSE, digits = 4)
