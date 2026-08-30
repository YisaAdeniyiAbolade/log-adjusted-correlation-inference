source("R/methods.R")
source("R/epa_data.R")

dir.create("results/epa", recursive = TRUE, showWarnings = FALSE)
dir.create("results/work", recursive = TRUE, showWarnings = FALSE)

data_2021 <- read_epa_data()$`2021`
platforms <- c("A10", "AIRPRO")
specs <- c("const", "u_add", "u_gain")

# Full-sample descriptive analysis.
full_rows <- list()
for (platform in platforms) {
  cfg <- platform_channels(platform, "2021")
  d <- complete_platform_data(data_2021, platform, "2021", TRUE)
  collocation <- sample_from_frame(d, seq_len(nrow(d)), platform, "2021", TRUE)
  target <- sample_from_frame(d, seq_len(nrow(d)), platform, "2021", FALSE)
  reference_rho <- cor(d[[cfg$xref]], d[[cfg$yref]])
  raw_rho <- cor(rowMeans(d[, cfg$w, drop = FALSE]), rowMeans(d[, cfg$z, drop = FALSE]))
  study_rho <- cor(rowMeans(d[, cfg$w_cal, drop = FALSE]), rowMeans(d[, cfg$z_cal, drop = FALSE]))

  for (spec in c("const", "u_add")) {
    est <- estimate_latent_correlation(collocation, target, spec)
    full_rows[[length(full_rows) + 1]] <- data.frame(
      platform = platform,
      platform_label = cfg$label,
      n = nrow(d),
      spec = spec,
      reference_rho = reference_rho,
      raw_sensor_rho = raw_rho,
      epa_study_calibrated_rho = study_rho,
      proxy_rho = est$proxy_rho,
      corrected_rho = est$rho,
      stabilized = est$stabilized,
      stringsAsFactors = FALSE
    )
  }
}
write.csv(do.call(rbind, full_rows), "results/epa/full_sample_results.csv", row.names = FALSE)

# Repeated held-out-burn evaluation. The reference values in the held-out set
# are used only to score the estimates, never to fit or compute the estimator.
set.seed(20260814)
reps <- 2000L
train_n <- 28L
base_burns <- complete_platform_data(data_2021, "A10", "2021")$burn
train_burns <- lapply(seq_len(reps), function(i) sample(base_burns, train_n, replace = FALSE))

heldout <- list()
k <- 1L
for (platform in platforms) {
  cfg <- platform_channels(platform, "2021")
  d <- complete_platform_data(data_2021, platform, "2021")
  if (!setequal(d$burn, base_burns)) stop("platforms do not share the same complete burn set")
  n <- nrow(d)

  for (rep in seq_len(reps)) {
    train <- which(d$burn %in% train_burns[[rep]])
    test <- setdiff(seq_len(n), train)
    collocation <- sample_from_frame(d, train, platform, "2021", TRUE)
    target <- sample_from_frame(d, test, platform, "2021", FALSE)
    benchmark <- cor(d[test, cfg$xref], d[test, cfg$yref])
    raw_rho <- cor(rowMeans(target$w), rowMeans(target$z))

    for (spec in specs) {
      rec <- data.frame(
        platform = platform,
        platform_label = cfg$label,
        rep = rep,
        spec = spec,
        n_collocation = length(train),
        n_heldout = length(test),
        target = benchmark,
        raw = raw_rho,
        proxy = NA_real_,
        corrected = NA_real_,
        stabilized = NA,
        feasible = FALSE,
        error = "",
        stringsAsFactors = FALSE
      )

      tryCatch({
        est <- estimate_latent_correlation(collocation, target, spec)
        rec$proxy <- est$proxy_rho
        rec$corrected <- est$rho
        rec$stabilized <- est$stabilized
        rec$feasible <- est$feasible
      }, error = function(e) {
        rec$error <<- substr(conditionMessage(e), 1, 160)
      })

      heldout[[k]] <- rec
      k <- k + 1L
    }
  }
}

heldout_results <- do.call(rbind, heldout)
write.csv(
  heldout_results,
  "results/work/heldout_burn_replications.csv",
  row.names = FALSE
)

summary_rows <- list()
k <- 1L
for (platform in platforms) {
  cfg <- platform_channels(platform, "2021")
  for (spec in specs) {
    g <- heldout_results[
      heldout_results$platform == platform & heldout_results$spec == spec,
      , drop = FALSE
    ]
    h <- g[complete.cases(g[, c("target", "raw", "proxy", "corrected")]), , drop = FALSE]

    for (method in c("raw", "proxy", "corrected")) {
      err <- h[[method]] - h$target
      ae <- abs(err)
      summary_rows[[k]] <- data.frame(
        platform = platform,
        platform_label = cfg$label,
        spec = spec,
        method = method,
        n_valid = nrow(h),
        bias = mean(err),
        mae = mean(ae),
        rmse = sqrt(mean(err^2)),
        median_ae = median(ae),
        p90_ae = as.numeric(quantile(ae, 0.90, names = FALSE)),
        stabilization_rate = mean(h$stabilized, na.rm = TRUE),
        feasible_rate = mean(h$feasible, na.rm = TRUE),
        failure_rate = mean(g$error != ""),
        fraction_improved = NA_real_,
        relative_mae_reduction = NA_real_,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }

    proxy_ae <- abs(h$proxy - h$target)
    corrected_ae <- abs(h$corrected - h$target)
    summary_rows[[k]] <- data.frame(
      platform = platform,
      platform_label = cfg$label,
      spec = spec,
      method = "corrected_vs_proxy",
      n_valid = nrow(h),
      bias = NA_real_,
      mae = NA_real_,
      rmse = NA_real_,
      median_ae = NA_real_,
      p90_ae = NA_real_,
      stabilization_rate = mean(h$stabilized, na.rm = TRUE),
      feasible_rate = mean(h$feasible, na.rm = TRUE),
      failure_rate = mean(g$error != ""),
      fraction_improved = mean(corrected_ae < proxy_ae),
      relative_mae_reduction = 1 - mean(corrected_ae) / mean(proxy_ae),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

write.csv(
  do.call(rbind, summary_rows),
  "results/epa/heldout_validation_summary.csv",
  row.names = FALSE
)
