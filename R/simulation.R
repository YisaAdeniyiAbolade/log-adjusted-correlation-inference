# Simulation mechanisms.

rmvnorm_base <- function(n, Sigma) {
  psd <- nearest_psd(Sigma)$matrix
  eig <- eigen(psd, symmetric = TRUE)
  A <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow = length(eig$values))
  matrix(rnorm(n * ncol(Sigma)), nrow = n) %*% t(A)
}

latent_parametric <- function(n, rho, family, params) {
  p <- params$latent_2021
  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)

  if (family == "normal") {
    q <- rmvnorm_base(n, Sigma)
  } else if (family == "t5") {
    df <- 5
    z <- rmvnorm_base(n, Sigma)
    scale <- sqrt(rchisq(n, df) / df)
    q <- z / scale * sqrt((df - 2) / df)
  } else if (family == "skewed") {
    slx <- 0.5
    sly <- 0.5
    denom <- sqrt((exp(slx^2) - 1) * (exp(sly^2) - 1))
    rlog <- log1p(rho * denom) / (slx * sly)
    if (abs(rlog) >= 0.999) stop("skewed correlation is not feasible")
    z <- rmvnorm_base(n, matrix(c(1, rlog, rlog, 1), 2, 2))
    ax <- exp(slx * z[, 1])
    ay <- exp(sly * z[, 2])
    qx <- (ax - exp(slx^2 / 2)) / sqrt((exp(slx^2) - 1) * exp(slx^2))
    qy <- (ay - exp(sly^2 / 2)) / sqrt((exp(sly^2) - 1) * exp(sly^2))
    q <- cbind(qx, qy)
  } else {
    stop("unknown latent distribution")
  }

  x <- p$mean_x + p$sd_x * q[, 1]
  y <- p$mean_y + p$sd_y * q[, 2]
  gamma <- params$rh_given_latent$gamma
  u <- gamma[1] + gamma[2] * q[, 1] + gamma[3] * q[, 2] +
    rnorm(n, 0, params$rh_given_latent$resid_sd)
  list(x = x, y = y, u = u)
}

latent_empirical <- function(n, params) {
  source_rows <- params$epa_2021_latent_rh
  out <- source_rows[sample(seq_len(nrow(source_rows)), n, replace = TRUE), , drop = FALSE]
  list(x = out[, 1], y = out[, 2], u = out[, 3])
}

error_covariance <- function(params, platform, structure = "full_epa", scale = 1,
                             transport_factor = 0) {
  S <- params$platforms[[platform]]$`2021`$raw_sensor_error_cov
  if (transport_factor > 0 && platform == "A10") {
    S19 <- params$platforms$A10$`2019`$raw_sensor_error_cov
    S <- (1 - transport_factor) * S + transport_factor * S19
  }
  if (structure == "diagonal") S <- diag(diag(S))
  if (structure == "no_cross_pollutant") {
    S[1:2, 3:4] <- 0
    S[3:4, 1:2] <- 0
  }
  if (!structure %in% c("full_epa", "diagonal", "no_cross_pollutant")) {
    stop("unknown error structure")
  }
  nearest_psd(S * scale^2)$matrix
}

generate_errors <- function(n, params, platform, structure, scale,
                            empirical = FALSE, transport_factor = 0) {
  if (empirical && structure == "full_epa" && transport_factor == 0) {
    E <- params$platforms[[platform]]$`2021`$raw_sensor_error_residuals
    out <- E[sample(seq_len(nrow(E)), n, replace = TRUE), , drop = FALSE]
    out <- sweep(out, 2, colMeans(E), "-")
    return(out * scale)
  }
  rmvnorm_base(n, error_covariance(params, platform, structure, scale, transport_factor))
}

make_readings <- function(x, y, u, errors, coefficients, spec = "u_add",
                          truth_variant = "u_add", quadratic_strength = 0.5) {
  latent <- cbind(x, x, y, y)
  readings <- matrix(NA_real_, length(x), 4)

  if (spec == "u_add") {
    for (j in 1:4) {
      readings[, j] <- coefficients[j, 1] + coefficients[j, 2] * u +
        coefficients[j, 3] * latent[, j] + errors[, j]
    }
  } else if (spec == "u_gain") {
    for (j in 1:4) {
      readings[, j] <- coefficients[j, 1] + coefficients[j, 2] * u +
        (coefficients[j, 3] + coefficients[j, 4] * u) * latent[, j] + errors[, j]
    }
  } else {
    stop("unknown systematic specification")
  }

  if (truth_variant == "quadratic_add") {
    uz <- (u - mean(u)) / max(sd(u), 1e-8)
    readings <- readings + quadratic_strength * outer(uz^2 - 1, apply(errors, 2, sd))
  }

  list(w = readings[, 1:2, drop = FALSE], z = readings[, 3:4, drop = FALSE])
}

generate_two_sample <- function(scenario, params) {
  value <- function(name, default = NULL) {
    x <- scenario[[name]][1]
    if (length(x) == 0 || is.na(x) || x == "") default else x
  }

  m <- as.integer(value("n_collocation"))
  n <- as.integer(value("n_target"))
  rho <- as.numeric(value("rho_true", params$latent_2021$rho))
  distribution <- as.character(value("latent_distribution", "normal"))
  truth_systematic <- as.character(value("truth_systematic", "u_add"))
  error_structure <- as.character(value("error_structure", "full_epa"))
  error_scale <- as.numeric(value("error_scale", 1))
  transport_shift <- as.character(value("transport_shift", "none"))
  platform <- as.character(value("platform", "A10"))

  if (distribution == "epa_empirical") {
    collocation_latent <- latent_empirical(m, params)
    target_latent <- latent_empirical(n, params)
    truth_rho <- params$latent_2021$rho
    empirical_error <- TRUE
  } else {
    collocation_latent <- latent_parametric(m, rho, distribution, params)
    target_latent <- latent_parametric(n, rho, distribution, params)
    truth_rho <- rho
    empirical_error <- FALSE
  }

  shift_factor <- switch(
    transport_shift,
    none = 0,
    moderate = 0.5,
    strong = 1,
    stop("unknown transport shift")
  )

  error_collocation <- generate_errors(
    m, params, platform, error_structure, error_scale, empirical_error, 0
  )
  error_target <- generate_errors(
    n, params, platform, error_structure, error_scale,
    empirical_error && shift_factor == 0, shift_factor
  )

  if (truth_systematic %in% c("u_add", "quadratic_add")) {
    coef_collocation <- raw_channel_coefficients(params, platform, "u_add", "2021")
    coef_target <- if (shift_factor == 0) {
      coef_collocation
    } else {
      transport_coefficients(params, platform, shift_factor)
    }
    read_collocation <- make_readings(
      collocation_latent$x, collocation_latent$y, collocation_latent$u,
      error_collocation, coef_collocation, "u_add", truth_systematic
    )
    read_target <- make_readings(
      target_latent$x, target_latent$y, target_latent$u,
      error_target, coef_target, "u_add", truth_systematic
    )
  } else if (truth_systematic == "u_gain") {
    coef_collocation <- raw_channel_coefficients(params, platform, "u_gain", "2021")
    read_collocation <- make_readings(
      collocation_latent$x, collocation_latent$y, collocation_latent$u,
      error_collocation, coef_collocation, "u_gain"
    )
    read_target <- make_readings(
      target_latent$x, target_latent$y, target_latent$u,
      error_target, coef_collocation, "u_gain"
    )
  } else {
    stop("unknown truth systematic specification")
  }

  collocation <- list(
    x = collocation_latent$x,
    y = collocation_latent$y,
    u = collocation_latent$u,
    w = read_collocation$w,
    z = read_collocation$z
  )
  target <- list(u = target_latent$u, w = read_target$w, z = read_target$z)
  truth <- list(rho = truth_rho, latent_sample_rho = cor(target_latent$x, target_latent$y))
  list(collocation = collocation, target = target, truth = truth)
}
