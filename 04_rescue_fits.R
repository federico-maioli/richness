library(tidyverse)
library(here)
library(sdmTMB)

# a quick pass over species that failed sanity() the first time - extra newton steps are
# much cheaper than a full refit and fix plain gradient/hessian issues, though not every
# cause (eg a genuinely unidentified spatial variance)
fit_log <- read_csv(here("results", "fit_log.csv"), show_col_types = FALSE)

failed <- fit_log |> filter(!converged, is.na(error))

cat("rescuing", nrow(failed), "species\n")

for (i in seq_len(nrow(failed))) {
  sp <- failed$species[i]
  group <- failed$group[i]
  path <- here("results", "fits", group, paste0(sp, ".rds"))

  fit <- readRDS(path)
  fit2 <- tryCatch(run_extra_optimization(fit, newton_loops = 3), error = \(e) e)

  if (inherits(fit2, "sdmTMB")) {
    ok <- sanity(fit2, gradient_thresh = 0.001)$all_ok
    if (ok) saveRDS(fit2, path)
  } else {
    ok <- FALSE
  }

  fit_log$converged[fit_log$species == sp] <- ok
  cat(group, sp, "rescued:", ok, "\n")
}

write_csv(fit_log, here("results", "fit_log.csv"))

cat("\ndone -", sum(fit_log$converged), "of", nrow(fit_log), "converged\n")
