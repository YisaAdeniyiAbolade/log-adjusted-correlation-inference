source("R/methods.R")
source("R/epa_data.R")
source("R/simulation.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: Rscript simulation/run_scenario.R SCENARIO_ID [REPS] [START_REP] [BOOTSTRAP_DRAWS]")
}

scenario_id <- args[1]
reps <- if (length(args) >= 2) as.integer(args[2]) else 1000L
start_rep <- if (length(args) >= 3) as.integer(args[3]) else 1L
bootstrap_B <- if (length(args) >= 4) as.integer(args[4]) else 499L
if (!is.finite(reps) || reps < 1) stop("REPS must be a positive integer")
if (!is.finite(start_rep) || start_rep < 1) stop("START_REP must be a positive integer")
if (!is.finite(bootstrap_B) || bootstrap_B < 19) stop("BOOTSTRAP_DRAWS must be at least 19")

scenarios <- read.csv("simulation/scenarios.csv", stringsAsFactors = FALSE)
scenario <- scenarios[scenarios$scenario_id == scenario_id, , drop = FALSE]
if (nrow(scenario) != 1) stop("scenario_id not found")
params <- build_epa_parameters()

ci_length <- function(ci) if (all(is.finite(ci))) ci[2] - ci[1] else NA_real_
raw_equal_rho <- function(target) cor(rowMeans(target$w), rowMeans(target$z))

one_rep <- function(rep_index, seed) {
  methods <- c("wald", "jel", "ajel", "bc_jel", "bc_ajel")
  true_default <- if (scenario$latent_distribution[1] == "epa_empirical") {
    params$latent_2021$rho
  } else {
    scenario$rho_true[1]
  }

  row <- list(
    rep = rep_index,
    seed = seed,
    rho_true = true_default,
    raw_equal_rho = NA_real_,
    proxy_rho = NA_real_,
    rho_hat = NA_real_,
    latent_sample_rho = NA_real_,
    raw_feasible = NA,
    stabilized = NA,
    point_error = NA_real_,
    proxy_error = NA_real_,
    point_failed = FALSE,
    inference_failed = TRUE,
    failure_message = "",
    jackknife_se = NA_real_,
    jel_lr_true = NA_real_,
    ajel_lr_true = NA_real_,
    bc_jel_critical = NA_real_,
    bc_ajel_critical = NA_real_
  )

  for (method in methods) {
    row[[paste0(method, "_low")]] <- NA_real_
    row[[paste0(method, "_high")]] <- NA_real_
    row[[paste0(method, "_length")]] <- NA_real_
    row[[paste0(method, "_cover")]] <- NA
  }

  set.seed(seed)
  stage <- tryCatch({
    generated <- generate_two_sample(scenario, params)
    collocation <- generated$collocation
    target <- generated$target
    spec <- scenario$fit_systematic[1]
    full <- estimate_latent_correlation(collocation, target, spec)
    list(
      collocation = collocation,
      target = target,
      truth = generated$truth,
      spec = spec,
      full = full
    )
  }, error = function(e) e)

  if (inherits(stage, "error")) {
    row$point_failed <- TRUE
    row$failure_message <- substr(conditionMessage(stage), 1, 240)
    return(as.data.frame(row, stringsAsFactors = FALSE))
  }

  collocation <- stage$collocation
  target <- stage$target
  truth <- stage$truth$rho
  full <- stage$full
  spec <- stage$spec

  row$rho_true <- truth
  row$raw_equal_rho <- raw_equal_rho(target)
  row$proxy_rho <- full$proxy_rho
  row$rho_hat <- full$rho
  row$latent_sample_rho <- stage$truth$latent_sample_rho
  row$raw_feasible <- full$feasible
  row$stabilized <- full$stabilized
  row$point_error <- full$rho - truth
  row$proxy_error <- full$proxy_rho - truth

  estimator <- function(v, m) estimate_latent_correlation(v, m, spec)$rho

  tryCatch({
    jk <- two_sample_jackknife(estimator, collocation, target)
    z <- variance_matched_pseudo(jk$theta, jk$collocation, jk$target)
    se <- jackknife_se(jk$collocation, jk$target)
    wald <- c(
      max(-1, jk$theta - qnorm(0.975) * se),
      min(1, jk$theta + qnorm(0.975) * se)
    )
    jel <- jel_ci(z, adjusted = FALSE)
    ajel <- jel_ci(z, adjusted = TRUE)

    bc_jel_critical <- influence_bootstrap_critical(
      jk$collocation, jk$target,
      B = bootstrap_B, seed = seed + 3000000L, adjusted = FALSE
    )
    bc_ajel_critical <- influence_bootstrap_critical(
      jk$collocation, jk$target,
      B = bootstrap_B, seed = seed + 6000000L, adjusted = TRUE
    )

    bc_jel <- if (is.finite(bc_jel_critical)) {
      jel_ci(z, adjusted = FALSE, critical_value = bc_jel_critical)
    } else {
      c(NA_real_, NA_real_)
    }
    bc_ajel <- if (is.finite(bc_ajel_critical)) {
      jel_ci(z, adjusted = TRUE, critical_value = bc_ajel_critical)
    } else {
      c(NA_real_, NA_real_)
    }

    intervals <- list(
      wald = wald,
      jel = jel,
      ajel = ajel,
      bc_jel = bc_jel,
      bc_ajel = bc_ajel
    )

    for (method in names(intervals)) {
      ci <- intervals[[method]]
      row[[paste0(method, "_low")]] <- ci[1]
      row[[paste0(method, "_high")]] <- ci[2]
      row[[paste0(method, "_length")]] <- ci_length(ci)
      row[[paste0(method, "_cover")]] <-
        all(is.finite(ci)) && ci[1] <= truth && truth <= ci[2]
    }

    row$jackknife_se <- se
    row$jel_lr_true <- jel_log_ratio(z, truth, FALSE)
    row$ajel_lr_true <- jel_log_ratio(z, truth, TRUE)
    row$bc_jel_critical <- bc_jel_critical
    row$bc_ajel_critical <- bc_ajel_critical
    row$inference_failed <- FALSE
  }, error = function(e) {
    row$failure_message <<- substr(conditionMessage(e), 1, 240)
  })

  as.data.frame(row, stringsAsFactors = FALSE)
}

scenario_index <- match(scenario_id, scenarios$scenario_id)
out <- vector("list", reps)
for (j in seq_len(reps)) {
  rep_index <- start_rep + j - 1L
  seed <- 20260814L + scenario_index * 100000L + rep_index
  out[[j]] <- one_rep(rep_index, seed)
}
out <- do.call(rbind, out)
for (name in names(scenario)) out[[name]] <- scenario[[name]][1]

dir.create("results/work/simulation", recursive = TRUE, showWarnings = FALSE)
write.csv(
  out,
  sprintf(
    "results/work/simulation/%s_%05d_%05d.csv",
    scenario_id, start_rep, start_rep + reps - 1L
  ),
  row.names = FALSE
)
