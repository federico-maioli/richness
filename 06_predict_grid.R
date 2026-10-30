library(tidyverse)
library(here)
library(sdmTMB)

set.seed(1)
n_sim <- 100

# only species whose fit actually converged - predict(nsim=...) draws from the joint
# precision matrix of the fixed and random effects, and that matrix is only meaningful
# if the model's own standard errors are. a species with an na or exploded se there
# would just inject noise into the richness draws rather than real uncertainty
fit_log <- read_csv(here("results", "fit_log.csv"), show_col_types = FALSE)
species_to_predict <- fit_log |> filter(converged) |> select(species, group)

cat("predicting", nrow(species_to_predict), "species onto the grid,", n_sim, "draws each\n")

# 01 Read the grid ----
# year has to be a factor matching the training labels for the (1 | year_factor) term -
# every grid year already exists in the training years, so no new levels show up
grid <- read_csv(here("data", "prediction_grid.csv"), show_col_types = FALSE) |>
  mutate(year_factor = factor(year)) |>
  as.data.frame()

# 02 Predict every species and accumulate richness ----
# a full species x node x year x draw array would be far too large to keep (278 species
# x 1,036,882 rows x 100 draws), so each species' simulated probabilities are added into
# a running richness total and then discarded - only the per species mean/sd are kept.
# other diversity metrics (shannon, etc) could be accumulated the same way inside this
# loop if needed later, without having to re-predict from scratch
richness_sims <- matrix(0, nrow = nrow(grid), ncol = n_sim)

dir.create(here("results", "predictions", "shelf"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "predictions", "deep"), recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(species_to_predict))) {
  sp <- species_to_predict$species[i]
  group <- species_to_predict$group[i]

  fit <- readRDS(here("results", "fits", group, paste0(sp, ".rds")))

  t <- system.time(
    sims <- predict(fit, newdata = grid, nsim = n_sim, type = "response")
  )

  richness_sims <- richness_sims + sims

  saveRDS(
    tibble(mean_prob = rowMeans(sims), sd_prob = apply(sims, 1, sd)),
    here("results", "predictions", group, paste0(sp, ".rds"))
  )

  cat(i, "/", nrow(species_to_predict), sp, round(t["elapsed"], 1), "sec\n")
}

# 03 Save the richness draws ----
richness_grid <- grid |>
  select(X, Y, lon, lat, year) |>
  mutate(
    richness_mean = rowMeans(richness_sims),
    richness_sd = apply(richness_sims, 1, sd),
    richness_q05 = apply(richness_sims, 1, \(x) quantile(x, 0.05)),
    richness_q95 = apply(richness_sims, 1, \(x) quantile(x, 0.95))
  )

write_csv(richness_grid, here("results", "richness_grid.csv"))
saveRDS(richness_sims, here("results", "richness_sims.rds"))

cat("\ndone - richness summary written to results/richness_grid.csv\n")
cat("raw draws (grid rows x", n_sim, "sims) saved to results/richness_sims.rds\n")
