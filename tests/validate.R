source("R/methods.R")
source("R/epa_data.R")
source("R/simulation.R")

# Data integrity checks for the tracked analysis extracts.
data <- read_epa_data()
stopifnot(nrow(data$`2019`) == 31, nrow(data$`2021`) == 40)
stopifnot(ncol(data$`2019`) == 8, ncol(data$`2021`) == 20)

params <- build_epa_parameters(data)
stopifnot(abs(params$latent_2021$rho - 0.8197466423) < 1e-8)

# Deterministic EPA full-sample checks under the primary RH-additive calibration.
check_platform <- function(platform, expected) {
  d <- complete_platform_data(data$`2021`, platform, "2021")
  collocation <- sample_from_frame(d, seq_len(nrow(d)), platform, "2021", TRUE)
  target <- sample_from_frame(d, seq_len(nrow(d)), platform, "2021", FALSE)
  est <- estimate_latent_correlation(collocation, target, "u_add")
  stopifnot(abs(est$rho - expected) < 1e-5)
  stopifnot(is.finite(est$proxy_rho), all(is.finite(est$noise_cov)))
}

check_platform("A10", 0.8061942953)
check_platform("AIRPRO", 0.7882141014)

# Exact recovery when channels obey an affine calibration with no random error.
u <- seq(-1, 1, length.out = 30)
t <- seq(0, 4, length.out = 30)
x <- 3 + 0.7 * cos(t) + 0.25 * t
y <- 2 + 0.6 * x + 0.4 * sin(1.7 * t)
w <- cbind(2 + 0.2 * u + 1.2 * x, -1 + 0.1 * u + 0.8 * x)
z <- cbind(1 - 0.3 * u + 1.1 * y, 0.5 + 0.15 * u + 0.9 * y)
collocation <- list(x = x, y = y, u = u, w = w, z = z)
target <- list(u = u, w = w, z = z)
est <- estimate_latent_correlation(collocation, target, "u_add")
stopifnot(abs(est$rho - cor(x, y)) < 1e-6)
stopifnot(max(abs(est$noise_cov)) < 1e-20)

# Positive-semidefinite safeguard and empirical-likelihood identities.
indefinite <- matrix(c(1, 2, 2, 1), 2, 2)
psd <- nearest_psd(indefinite)
stopifnot(psd$changed, min(eigen(psd$matrix, symmetric = TRUE)$values) > 0)

ztest <- c(-2, -1, 0, 1, 2) + 0.25
stopifnot(abs(jel_log_ratio(ztest, mean(ztest))) < 1e-10)
stopifnot(is.finite(jel_log_ratio(ztest, 0.25, adjusted = TRUE)))

# One generated two-sample scenario exercises the full point and jackknife path.
scenario <- read.csv("simulation/scenarios.csv", stringsAsFactors = FALSE)[1, , drop = FALSE]
set.seed(20260814)
generated <- generate_two_sample(scenario, params)
point <- estimate_latent_correlation(
  generated$collocation,
  generated$target,
  scenario$fit_systematic[1]
)
stopifnot(is.finite(point$rho))

estimator <- function(v, m) {
  estimate_latent_correlation(v, m, scenario$fit_systematic[1])$rho
}
jk <- two_sample_jackknife(estimator, generated$collocation, generated$target)
zpool <- variance_matched_pseudo(jk$theta, jk$collocation, jk$target)
stopifnot(is.finite(jackknife_se(jk$collocation, jk$target)))
stopifnot(abs(mean(zpool) - jk$theta) < 1e-10)

cat("All validation checks passed\n")

# Reference-free log-adjusted checks.
u_lar <- seq(-1, 1, length.out = 40)
t_lar <- seq(0, 2 * pi, length.out = 40)
x_lar <- 0.7 * cos(t_lar) + 0.2 * t_lar
y_lar <- 0.55 * x_lar + 0.35 * sin(1.5 * t_lar)
w_lar <- x_lar + 0.6 * u_lar
z_lar <- y_lar - 0.4 * u_lar
lar_fit <- lar_correlation(w_lar, z_lar, u_lar)
melar_zero <- melar_correlation(w_lar, z_lar, u_lar, matrix(0, 2, 2))
stopifnot(is.finite(lar_fit$rho), abs(lar_fit$rho - melar_zero$rho) < 1e-12)
lar_jk <- lar_melar_jackknife(w_lar, z_lar, u_lar, method = "lar")
stopifnot(is.finite(lar_jk$theta), all(is.finite(lar_jk$pseudo)))
