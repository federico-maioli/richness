library(tidyverse)
library(here)
library(sdmTMB)
library(pROC)

# in sample auc and true skill statistic for every fitted species. predict() carries the
# original presence column through, so no need to reload the training data separately
fit_log <- read_csv(here("results", "fit_log.csv"), show_col_types = FALSE)

skill <- tibble(
  species = character(), group = character(), converged = logical(),
  auc = numeric(),
  tss_fixed = numeric(), true_pos_rate_fixed = numeric(), true_neg_rate_fixed = numeric(),
  tss_optimal = numeric(), threshold_optimal = numeric()
)

for (i in seq_len(nrow(fit_log))) {
  sp <- fit_log$species[i]
  group <- fit_log$group[i]

  path <- here("results", "fits", group, paste0(sp, ".rds"))
  if (!file.exists(path)) next

  fit <- readRDS(path)
  pred <- predict(fit)

  actual <- pred[[sp]]
  prob <- plogis(pred$est)

  roc <- pROC::roc(actual, prob, quiet = TRUE)
  auc <- as.numeric(pROC::auc(roc))

  # true skill statistic at the standard 0.5 threshold - a fixed threshold like this one
  # can read as 0 for a rare species whose probabilities never cross 0.5, even when the
  # model ranks presences well, so it is kept alongside the threshold free version below
  predicted_class <- ifelse(prob < 0.5, 0, 1)
  conmat <- table(factor(predicted_class, levels = c(0, 1)), factor(actual, levels = c(0, 1)))

  true_neg <- conmat[1, 1]
  false_neg <- conmat[1, 2]
  false_pos <- conmat[2, 1]
  true_pos <- conmat[2, 2]

  true_pos_rate_fixed <- true_pos / (true_pos + false_neg)
  true_neg_rate_fixed <- true_neg / (true_neg + false_pos)
  tss_fixed <- true_pos_rate_fixed + true_neg_rate_fixed - 1

  # tss at its own best threshold instead - the youden index, which is what "best"
  # finds, is sensitivity + specificity - 1, the same quantity as tss
  best <- pROC::coords(roc, "best", best.method = "youden", transpose = FALSE)
  tss_optimal <- best$sensitivity[1] + best$specificity[1] - 1
  threshold_optimal <- best$threshold[1]

  skill <- bind_rows(skill, tibble(
    species = sp, group = group, converged = fit_log$converged[i],
    auc = auc,
    tss_fixed = tss_fixed, true_pos_rate_fixed = true_pos_rate_fixed, true_neg_rate_fixed = true_neg_rate_fixed,
    tss_optimal = tss_optimal, threshold_optimal = threshold_optimal
  ))

  cat(sp, "auc:", round(auc, 3), "tss fixed:", round(tss_fixed, 3), "tss optimal:", round(tss_optimal, 3), "\n")
}

write_csv(skill, here("results", "predictive_skill.csv"))

cat("\nsummary, converged species only:\n")
print(skill |> filter(converged) |> summarise(
  median_auc = median(auc), median_tss_fixed = median(tss_fixed), median_tss_optimal = median(tss_optimal),
  min_auc = min(auc), min_tss_optimal = min(tss_optimal)
))

cat("\nsummary, all species:\n")
print(skill |> summarise(
  median_auc = median(auc), median_tss_fixed = median(tss_fixed), median_tss_optimal = median(tss_optimal),
  min_auc = min(auc), min_tss_optimal = min(tss_optimal)
))
