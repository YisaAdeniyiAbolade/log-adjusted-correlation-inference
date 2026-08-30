args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 1000L
bootstrap_draws <- if (length(args) >= 2) as.integer(args[2]) else 499L
if (!is.finite(reps) || reps < 1) stop("REPS must be a positive integer")
if (!is.finite(bootstrap_draws) || bootstrap_draws < 19) {
  stop("BOOTSTRAP_DRAWS must be at least 19")
}

scenarios <- read.csv("simulation/scenarios.csv", stringsAsFactors = FALSE)
dir.create("results/work/simulation", recursive = TRUE, showWarnings = FALSE)
old <- list.files("results/work/simulation", pattern = "\\.csv$", full.names = TRUE)
if (length(old)) file.remove(old)

for (id in scenarios$scenario_id) {
  message("Running scenario ", id)
  status <- system2(
    "Rscript",
    c("simulation/run_scenario.R", id, reps, 1L, bootstrap_draws)
  )
  if (!identical(status, 0L)) stop(sprintf("scenario %s failed", id))
}

source("simulation/summarize_results.R")
