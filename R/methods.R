# Core methods for log-adjusted and calibration-assisted correlation inference.
# The R implementation uses base R only.

epanechnikov_kernel <- function(t) {
  t <- as.numeric(t)
  out <- 0.75 * (1 - t^2)
  out[abs(t) > 1] <- 0
  out
}

log_mean_exp <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || any(!is.finite(x))) stop("x must contain finite values")
  a <- max(x)
  a + log(mean(exp(x - a)))
}

lar_adjust <- function(values, u, bandwidth = NULL) {
  values <- as.numeric(values)
  u <- as.numeric(u)
  n <- length(values)
  if (length(u) != n || n < 3 || any(!is.finite(values)) || any(!is.finite(u))) {
    stop("values and u must be finite vectors of equal length >= 3")
  }
  if (is.null(bandwidth)) bandwidth <- max(sd(u) * n^(-1/3), 1e-8)
  if (!is.finite(bandwidth) || bandwidth <= 0) stop("bandwidth must be positive")
  conditional_log_mean <- numeric(n)
  for (i in seq_len(n)) {
    weights <- epanechnikov_kernel((u - u[i]) / bandwidth)
    den <- sum(weights)
    if (!is.finite(den) || den <= 0) stop("kernel denominator is zero")
    a <- max(values)
    conditional_log_mean[i] <- a + log(sum(weights * exp(values - a)) / den)
  }
  values - conditional_log_mean + log_mean_exp(values)
}

corrected_log_correlation <- function(x, y, error_cov = matrix(0, 2, 2), clip = TRUE) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(x) != length(y) || length(x) < 3) stop("x and y must have equal length >= 3")
  error_cov <- as.matrix(error_cov)
  if (!all(dim(error_cov) == c(2, 2)) || any(!is.finite(error_cov))) {
    stop("error_cov must be a finite 2 x 2 matrix")
  }
  latent_cov <- cov(cbind(x, y)) - error_cov
  latent_cov <- (latent_cov + t(latent_cov)) / 2
  if (latent_cov[1, 1] <= 1e-12 || latent_cov[2, 2] <= 1e-12) return(NA_real_)
  rho <- latent_cov[1, 2] / sqrt(latent_cov[1, 1] * latent_cov[2, 2])
  if (clip) rho <- max(-0.999999, min(0.999999, rho))
  rho
}

lar_correlation <- function(w, z, u, bandwidth = NULL) {
  x_adjusted <- lar_adjust(w, u, bandwidth)
  y_adjusted <- lar_adjust(z, u, bandwidth)
  list(
    rho = corrected_log_correlation(x_adjusted, y_adjusted),
    x_adjusted = x_adjusted,
    y_adjusted = y_adjusted
  )
}

melar_correlation <- function(w, z, u, error_cov, bandwidth = NULL) {
  x_adjusted <- lar_adjust(w, u, bandwidth)
  y_adjusted <- lar_adjust(z, u, bandwidth)
  list(
    rho = corrected_log_correlation(x_adjusted, y_adjusted, error_cov),
    x_adjusted = x_adjusted,
    y_adjusted = y_adjusted
  )
}

lar_melar_jackknife <- function(w, z, u, error_cov = matrix(0, 2, 2),
                                method = c("melar", "lar"), bandwidth = NULL) {
  method <- match.arg(method)
  w <- as.numeric(w); z <- as.numeric(z); u <- as.numeric(u)
  n <- length(w)
  if (length(z) != n || length(u) != n || n < 4) stop("w, z and u must have equal length >= 4")
  if (is.null(bandwidth)) bandwidth <- max(sd(u) * n^(-1/3), 1e-8)
  estimate <- function(wi, zi, ui) {
    if (method == "lar") lar_correlation(wi, zi, ui, bandwidth)$rho
    else melar_correlation(wi, zi, ui, error_cov, bandwidth)$rho
  }
  theta <- estimate(w, z, u)
  leave <- numeric(n)
  for (i in seq_len(n)) {
    keep <- seq_len(n) != i
    leave[i] <- estimate(w[keep], z[keep], u[keep])
  }
  if (!is.finite(theta) || any(!is.finite(leave))) {
    return(list(theta = theta, pseudo = rep(NA_real_, n), bandwidth = bandwidth))
  }
  list(theta = theta, pseudo = n * theta - (n - 1) * leave, bandwidth = bandwidth)
}

normalize_weights <- function(weights, k, label = "weights") {
  if (is.null(weights)) return(rep(1 / k, k))
  weights <- as.numeric(weights)
  if (length(weights) != k || any(!is.finite(weights))) {
    stop(sprintf("%s must contain %d finite values", label, k))
  }
  total <- sum(weights)
  if (!is.finite(total) || abs(total) < 1e-12) {
    stop(sprintf("%s must have a nonzero sum", label))
  }
  weights / total
}

fit_calibration <- function(x, recorded, u, spec = c("u_add", "const", "u_gain")) {
  spec <- match.arg(spec)
  x <- as.numeric(x)
  recorded <- as.numeric(recorded)
  u <- as.numeric(u)
  n <- length(x)
  if (length(recorded) != n || length(u) != n || n < 3) {
    stop("calibration inputs must have equal length and at least three observations")
  }
  if (any(!is.finite(x)) || any(!is.finite(recorded)) || any(!is.finite(u))) {
    stop("calibration inputs must be finite")
  }
  u0 <- mean(u)
  us <- sd(u)
  if (!is.finite(us) || us <= 0) stop("calibration covariate has no variation")
  h <- (u - u0) / us
  design <- switch(
    spec,
    const = cbind(1, x),
    u_add = cbind(1, h, x),
    u_gain = cbind(1, h, x, x * h)
  )
  fit <- lm.fit(design, recorded)
  if (fit$rank < ncol(design) || any(!is.finite(fit$coefficients))) {
    stop("calibration fit is rank deficient")
  }
  list(beta = fit$coefficients, u0 = u0, us = us, spec = spec)
}

calibration_functions <- function(fit, u) {
  u <- as.numeric(u)
  h <- (u - fit$u0) / fit$us
  b <- fit$beta
  if (fit$spec == "const") {
    return(list(delta = rep(b[1], length(u)), phi = rep(b[2], length(u))))
  }
  if (fit$spec == "u_add") {
    return(list(delta = b[1] + b[2] * h, phi = rep(b[3], length(u))))
  }
  list(delta = b[1] + b[2] * h, phi = b[3] + b[4] * h)
}

inverse_calibrate <- function(recorded, u, fit, gain_tol = 1e-8) {
  functions <- calibration_functions(fit, u)
  if (any(!is.finite(functions$phi)) || any(abs(functions$phi) < gain_tol)) {
    stop("fitted multiplicative gain is too close to zero")
  }
  (as.numeric(recorded) - functions$delta) / functions$phi
}

fit_calibration_bundle <- function(x, y, u, w, z, spec = c("u_add", "const", "u_gain")) {
  spec <- match.arg(spec)
  w <- as.matrix(w)
  z <- as.matrix(z)
  if (nrow(w) != length(x) || nrow(z) != length(y) || length(x) != length(y) || length(u) != length(x)) {
    stop("collocation variables and channel matrices have incompatible dimensions")
  }
  fits_x <- lapply(seq_len(ncol(w)), function(j) fit_calibration(x, w[, j], u, spec))
  fits_y <- lapply(seq_len(ncol(z)), function(j) fit_calibration(y, z[, j], u, spec))
  list(fits_x = fits_x, fits_y = fits_y, spec = spec)
}

apply_calibration_bundle <- function(bundle, u, w, z) {
  w <- as.matrix(w)
  z <- as.matrix(z)
  if (nrow(w) != length(u) || nrow(z) != length(u)) {
    stop("target covariate and channel matrices have incompatible dimensions")
  }
  if (ncol(w) != length(bundle$fits_x) || ncol(z) != length(bundle$fits_y)) {
    stop("channel counts do not match the fitted calibration bundle")
  }
  px <- sapply(seq_len(ncol(w)), function(j) inverse_calibrate(w[, j], u, bundle$fits_x[[j]]))
  py <- sapply(seq_len(ncol(z)), function(j) inverse_calibrate(z[, j], u, bundle$fits_y[[j]]))
  if (is.null(dim(px))) px <- matrix(px, ncol = 1)
  if (is.null(dim(py))) py <- matrix(py, ncol = 1)
  list(x = px, y = py)
}

aggregate_channels <- function(channel_matrix, weights = NULL) {
  channel_matrix <- as.matrix(channel_matrix)
  weights <- normalize_weights(weights, ncol(channel_matrix), "channel weights")
  drop(channel_matrix %*% weights)
}

aggregate_error_covariance <- function(error_matrix, rx, ry, wx = NULL, wy = NULL) {
  error_matrix <- as.matrix(error_matrix)
  if (ncol(error_matrix) != rx + ry) stop("error matrix has an unexpected number of columns")
  wx <- normalize_weights(wx, rx, "X-channel weights")
  wy <- normalize_weights(wy, ry, "Y-channel weights")
  A <- matrix(0, 2, rx + ry)
  A[1, seq_len(rx)] <- wx
  A[2, rx + seq_len(ry)] <- wy
  A %*% cov(error_matrix) %*% t(A)
}

estimate_collocation_noise <- function(collocation, bundle, wx = NULL, wy = NULL) {
  proxy <- apply_calibration_bundle(bundle, collocation$u, collocation$w, collocation$z)
  ex <- sweep(proxy$x, 1, collocation$x, "-")
  ey <- sweep(proxy$y, 1, collocation$y, "-")
  error_matrix <- cbind(ex, ey)
  aggregate_error_covariance(error_matrix, ncol(ex), ncol(ey), wx, wy)
}

nearest_psd <- function(S, eps = 1e-10) {
  S <- as.matrix(S)
  if (nrow(S) != ncol(S) || nrow(S) < 1 || any(!is.finite(S))) {
    stop("S must be a finite square matrix")
  }
  S <- (S + t(S)) / 2
  eig <- eigen(S, symmetric = TRUE)
  scale <- max(sum(abs(diag(S))), 1)
  floor <- eps * scale
  changed <- min(eig$values) < floor
  values <- pmax(eig$values, floor)
  out <- eig$vectors %*% diag(values, nrow = length(values)) %*% t(eig$vectors)
  out <- (out + t(out)) / 2
  list(matrix = out, changed = changed, eigenvalues_raw = eig$values, floor = floor)
}

estimate_latent_correlation <- function(collocation, target,
                                        spec = c("u_add", "const", "u_gain"),
                                        wx = NULL, wy = NULL,
                                        stabilize = TRUE, eps = 1e-10) {
  spec <- match.arg(spec)
  m <- length(collocation$u)
  n <- length(target$u)
  if (m < 3 || n < 2) stop("insufficient collocation or target observations")
  required_collocation <- c("x", "y", "u", "w", "z")
  required_target <- c("u", "w", "z")
  if (!all(required_collocation %in% names(collocation))) stop("collocation sample is incomplete")
  if (!all(required_target %in% names(target))) stop("target sample is incomplete")

  bundle <- fit_calibration_bundle(
    collocation$x, collocation$y, collocation$u,
    collocation$w, collocation$z, spec
  )
  noise_cov <- estimate_collocation_noise(collocation, bundle, wx, wy)
  proxy_channels <- apply_calibration_bundle(bundle, target$u, target$w, target$z)
  xp <- aggregate_channels(proxy_channels$x, wx)
  yp <- aggregate_channels(proxy_channels$y, wy)
  observed_cov <- cov(cbind(xp, yp))
  latent_raw <- (observed_cov - noise_cov + t(observed_cov - noise_cov)) / 2

  feasible <- is.finite(latent_raw[1, 1]) && is.finite(latent_raw[2, 2]) &&
    latent_raw[1, 1] > 0 && latent_raw[2, 2] > 0 &&
    is.finite(latent_raw[1, 2]) &&
    abs(latent_raw[1, 2]) < sqrt(latent_raw[1, 1] * latent_raw[2, 2])

  stabilized <- FALSE
  latent_cov <- latent_raw
  if (stabilize) {
    psd <- nearest_psd(latent_raw, eps)
    latent_cov <- psd$matrix
    stabilized <- psd$changed
  } else if (!feasible) {
    return(list(
      rho = NA_real_, proxy_rho = cor(xp, yp), latent_cov = latent_raw,
      latent_cov_raw = latent_raw, observed_cov = observed_cov, noise_cov = noise_cov,
      stabilized = FALSE, feasible = FALSE, calibration = bundle
    ))
  }

  rho <- latent_cov[1, 2] / sqrt(latent_cov[1, 1] * latent_cov[2, 2])
  rho <- max(-1, min(1, rho))
  list(
    rho = rho,
    proxy_rho = cor(xp, yp),
    latent_cov = latent_cov,
    latent_cov_raw = latent_raw,
    observed_cov = observed_cov,
    noise_cov = noise_cov,
    stabilized = stabilized,
    feasible = feasible,
    calibration = bundle
  )
}

subset_sample <- function(sample, keep) {
  out <- sample
  n <- length(sample$u)
  for (nm in names(sample)) {
    value <- sample[[nm]]
    if (is.matrix(value) && nrow(value) == n) out[[nm]] <- value[keep, , drop = FALSE]
    if (!is.matrix(value) && length(value) == n) out[[nm]] <- value[keep]
  }
  out
}

two_sample_jackknife <- function(estimator, collocation, target) {
  m <- length(collocation$u)
  n <- length(target$u)
  theta <- estimator(collocation, target)
  if (!is.finite(theta)) stop("full estimator is not finite")

  pv_target <- numeric(n)
  for (j in seq_len(n)) {
    keep <- seq_len(n) != j
    theta_minus <- estimator(collocation, subset_sample(target, keep))
    if (!is.finite(theta_minus)) stop("non-finite target delete-one estimate")
    pv_target[j] <- n * theta - (n - 1) * theta_minus
  }

  pv_collocation <- numeric(m)
  for (i in seq_len(m)) {
    keep <- seq_len(m) != i
    theta_minus <- estimator(subset_sample(collocation, keep), target)
    if (!is.finite(theta_minus)) stop("non-finite collocation delete-one estimate")
    pv_collocation[i] <- m * theta - (m - 1) * theta_minus
  }

  list(theta = theta, collocation = pv_collocation, target = pv_target)
}

variance_matched_pseudo <- function(theta, pv_collocation, pv_target) {
  m <- length(pv_collocation)
  n <- length(pv_target)
  N <- m + n
  zv <- theta + (N / m) * (pv_collocation - mean(pv_collocation))
  zm <- theta + (N / n) * (pv_target - mean(pv_target))
  z <- c(zv, zm)
  # Restore the pooled mean to theta up to numerical precision.
  z + theta - mean(z)
}

jackknife_se <- function(pv_collocation, pv_target) {
  sqrt(var(pv_collocation) / length(pv_collocation) + var(pv_target) / length(pv_target))
}

el_lambda <- function(g) {
  g <- as.numeric(g)
  if (all(abs(g) < 1e-12)) return(0)
  gmin <- min(g)
  gmax <- max(g)
  if (gmin >= 0 || gmax <= 0) return(NA_real_)
  lo <- -1 / gmax + 1e-11
  hi <- -1 / gmin - 1e-11
  score <- function(lambda) sum(g / (1 + lambda * g))
  tryCatch(uniroot(score, c(lo, hi), maxiter = 200)$root, error = function(e) NA_real_)
}

jel_log_ratio <- function(z, theta, adjusted = FALSE, a = NULL) {
  g <- as.numeric(z) - theta
  if (adjusted) {
    N <- length(g)
    if (is.null(a)) a <- log(max(N, 2)) / 2
    g <- c(g, -a * mean(g))
  }
  lambda <- el_lambda(g)
  if (!is.finite(lambda)) return(Inf)
  den <- 1 + lambda * g
  if (any(den <= 0)) return(Inf)
  2 * sum(log(den))
}

jel_ci <- function(z, level = 0.95, adjusted = FALSE, critical_value = NULL,
                   domain = c(-0.999, 0.999), grid_size = 401) {
  crit <- if (is.null(critical_value)) qchisq(level, 1) else critical_value
  grid <- seq(domain[1], domain[2], length.out = grid_size)
  values <- vapply(grid, function(theta) jel_log_ratio(z, theta, adjusted), numeric(1))
  ok <- is.finite(values) & values <= crit
  if (!any(ok)) return(c(NA_real_, NA_real_))

  ids <- which(ok)
  li <- min(ids)
  ri <- max(ids)
  left <- grid[li]
  right <- grid[ri]
  root_fun <- function(theta) jel_log_ratio(z, theta, adjusted) - crit

  if (li > 1 && is.finite(values[li - 1])) {
    left <- tryCatch(uniroot(root_fun, c(grid[li - 1], grid[li]))$root, error = function(e) left)
  }
  if (ri < length(grid) && is.finite(values[ri + 1])) {
    right <- tryCatch(uniroot(root_fun, c(grid[ri], grid[ri + 1]))$root, error = function(e) right)
  }
  c(left, right)
}

influence_bootstrap_critical <- function(pv_collocation, pv_target, level = 0.95,
                                         B = 199, seed = 20260814,
                                         adjusted = TRUE) {
  if (length(pv_collocation) < 3 || length(pv_target) < 3) {
    stop("insufficient pseudo-values")
  }
  if (!is.finite(B) || B < 19) stop("B must be at least 19")
  set.seed(seed)
  m <- length(pv_collocation)
  n <- length(pv_target)
  N <- m + n
  cv <- pv_collocation - mean(pv_collocation)
  ct <- pv_target - mean(pv_target)
  lr <- numeric(B)

  for (b in seq_len(B)) {
    vb <- sample(cv, m, replace = TRUE)
    tb <- sample(ct, n, replace = TRUE)
    delta <- mean(vb) + mean(tb)
    zstar <- c(
      delta + (N / m) * (vb - mean(vb)),
      delta + (N / n) * (tb - mean(tb))
    )
    lr[b] <- jel_log_ratio(zstar, 0, adjusted)
  }

  finite <- lr[is.finite(lr)]
  if (length(finite) < 0.8 * B) return(NA_real_)
  as.numeric(quantile(finite, level, names = FALSE, type = 7))
}
