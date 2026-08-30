dir.create("results/simulation", recursive = TRUE, showWarnings = FALSE)
files <- list.files("results/work/simulation", pattern = "\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("no simulation result files found")

all_results <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
scenarios <- read.csv("simulation/scenarios.csv", stringsAsFactors = FALSE)
methods <- c("wald", "jel", "ajel", "bc_jel", "bc_ajel")
summary_rows <- list()

for (i in seq_len(nrow(scenarios))) {
  id <- scenarios$scenario_id[i]
  g <- all_results[all_results$scenario_id == id, , drop = FALSE]
  if (nrow(g) == 0) next
  q <- g[!as.logical(g$inference_failed) & !as.logical(g$point_failed), , drop = FALSE]

  row <- as.list(scenarios[i, ])
  row$reps <- nrow(g)
  row$point_mean <- mean(g$rho_hat, na.rm = TRUE)
  row$point_bias <- mean(g$point_error, na.rm = TRUE)
  row$point_rmse <- sqrt(mean(g$point_error^2, na.rm = TRUE))
  row$proxy_mean <- mean(g$proxy_rho, na.rm = TRUE)
  row$proxy_bias <- mean(g$proxy_error, na.rm = TRUE)
  row$proxy_rmse <- sqrt(mean(g$proxy_error^2, na.rm = TRUE))
  row$stabilization_rate <- mean(as.logical(g$stabilized), na.rm = TRUE)
  row$raw_infeasible_rate <- 1 - mean(as.logical(g$raw_feasible), na.rm = TRUE)
  row$point_failure_rate <- mean(as.logical(g$point_failed), na.rm = TRUE)
  row$inference_failure_rate <- mean(as.logical(g$inference_failed), na.rm = TRUE)

  for (method in methods) {
    length_col <- paste0(method, "_length")
    cover_col <- paste0(method, "_cover")
    finite <- is.finite(q[[length_col]])
    row[[paste0(method, "_coverage")]] <- if (any(finite)) {
      mean(as.logical(q[[cover_col]][finite]), na.rm = TRUE)
    } else {
      NA_real_
    }
    row[[paste0(method, "_avg_length")]] <- if (any(finite)) {
      mean(q[[length_col]][finite], na.rm = TRUE)
    } else {
      NA_real_
    }
    row[[paste0(method, "_finite_rate")]] <- if (nrow(q) > 0) mean(finite) else NA_real_
  }

  summary_rows[[length(summary_rows) + 1]] <- as.data.frame(row, stringsAsFactors = FALSE)
}

summary <- do.call(rbind, summary_rows)
write.csv(
  summary[summary$scenario_family %in% c("EPA_primary", "EPA_replication"), ],
  "results/simulation/primary_summary.csv",
  row.names = FALSE
)
write.csv(
  summary[!summary$scenario_family %in% c("EPA_primary", "EPA_replication"), ],
  "results/simulation/robustness_summary.csv",
  row.names = FALSE
)
