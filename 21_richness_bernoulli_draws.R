library(tidyverse)
library(here)

# richness_sims.rds (from 06_predict_grid.R) sums each species' simulated PROBABILITY
# across species - that captures model/parameter uncertainty (the draws come from the
# joint precision matrix) but still treats richness as a continuous "expected" count.
# this instead treats each species' predicted probability as a bernoulli distribution -
# present or absent, nothing in between - draws from it 100 times per species, and
# sums across species per draw. that adds the observation level uncertainty a real
# survey would show (would this exact species actually have been caught here), on top
# of the point estimate probability already fitted. it does not re-add the earlier
# parameter uncertainty, since that would mean re-predicting all 221 species with
# nsim draws again (another ~50 minute run) rather than reusing the saved point
# estimates - ask if that combined version is wanted
set.seed(1)
n_draws <- 100

fit_log <- read_csv(here("results", "fit_log.csv"), show_col_types = FALSE)
species_to_use <- fit_log |> filter(converged) |> select(species, group)

cat("combining", nrow(species_to_use), "species,", n_draws, "draws each\n")

n_rows <- nrow(read_csv(here("results", "richness_grid.csv"), show_col_types = FALSE, col_select = "X"))

richness_draws <- matrix(0L, nrow = n_rows, ncol = n_draws)

for (i in seq_len(nrow(species_to_use))) {
  sp <- species_to_use$species[i]
  group <- species_to_use$group[i]

  mean_prob <- readRDS(here("results", "predictions", group, paste0(sp, ".rds")))$mean_prob

  presence <- matrix(
    rbinom(n_rows * n_draws, size = 1, prob = mean_prob),
    nrow = n_rows, ncol = n_draws
  )

  richness_draws <- richness_draws + presence

  if (i %% 25 == 0) cat(i, "/", nrow(species_to_use), "\n")
}

cat("\nbernoulli richness draws done\n")

# 02 Save a summary, same shape as richness_grid.csv, plus the raw draws ----
grid_xy <- read_csv(here("results", "richness_grid.csv"), show_col_types = FALSE) |>
  select(X, Y, lon, lat, year)

richness_bernoulli <- grid_xy |>
  mutate(
    richness_mean = rowMeans(richness_draws),
    richness_sd = apply(richness_draws, 1, sd),
    richness_q05 = apply(richness_draws, 1, \(x) quantile(x, 0.05)),
    richness_q95 = apply(richness_draws, 1, \(x) quantile(x, 0.95))
  )

write_csv(richness_bernoulli, here("results", "richness_bernoulli.csv"))
saveRDS(richness_draws, here("results", "richness_bernoulli_draws.rds"))

# 03 Compare the spread against the probability sum version already on disk ----
richness_prob <- read_csv(here("results", "richness_grid.csv"), show_col_types = FALSE)

cat("\nprobability sum richness (06_predict_grid.R) - sd summary:\n")
print(summary(richness_prob$richness_sd))

cat("\nbernoulli richness - sd summary:\n")
print(summary(richness_bernoulli$richness_sd))

cat(
  "\nmean richness, probability sum version:", round(mean(richness_prob$richness_mean), 2),
  "| bernoulli version:", round(mean(richness_bernoulli$richness_mean), 2), "\n"
)

cat("\nwritten to results/richness_bernoulli.csv and results/richness_bernoulli_draws.rds\n")
