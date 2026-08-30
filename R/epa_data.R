read_epa_data <- function() {
  list(
    `2019` = read.csv("data/epa_2019_burn_integrated.csv", check.names = FALSE),
    `2021` = read.csv("data/epa_2021_burn_integrated.csv", check.names = FALSE)
  )
}

platform_channels <- function(platform, year) {
  if (platform == "A10" && year == "2021") {
    return(list(
      label = "Kunak Model A10",
      xref = "Reference_CO", yref = "Reference_NO2",
      w = c("kunak_A14_1_CO", "kunak_A14_2_CO"),
      z = c("kunak_A14_1_NO2", "kunak_A14_2_NO2"),
      w_cal = c("kunak_A14_1_CO_lin", "kunak_A14_2_CO_lin"),
      z_cal = c("kunak_A14_1_NO2_lin", "kunak_A14_2_NO2_lin")
    ))
  }
  if (platform == "A10" && year == "2019") {
    return(list(
      label = "Kunak Model A10",
      xref = "Reference_CO", yref = "Reference_NO2",
      w = c("kunak_1_CO_ppm", "kunak_2_CO_ppm"),
      z = c("kunak_1_NO2_ppb", "kunak_2_NO2_ppb")
    ))
  }
  if (platform == "AIRPRO" && year == "2021") {
    return(list(
      label = "Kunak Air Pro",
      xref = "Reference_CO", yref = "Reference_NO2",
      w = c("kunak_A3_1_CO", "kunak_A3_2_CO"),
      z = c("kunak_A3_1_NO2", "kunak_A3_2_NO2"),
      w_cal = c("kunak_A3_1_CO_linftz", "kunak_A3_2_CO_linftz"),
      z_cal = c("kunak_A3_1_NO2_quad", "kunak_A3_2_NO2_quad")
    ))
  }
  stop("platform/campaign combination is not available")
}

complete_platform_data <- function(data, platform, year, include_study_calibration = FALSE) {
  cfg <- platform_channels(platform, year)
  cols <- c("burn", "chamber_RH", cfg$xref, cfg$yref, cfg$w, cfg$z)
  if (include_study_calibration && !is.null(cfg$w_cal)) cols <- c(cols, cfg$w_cal, cfg$z_cal)
  out <- data[, cols, drop = FALSE]
  out[complete.cases(out), , drop = FALSE]
}

sample_from_frame <- function(data, idx, platform, year, references = TRUE) {
  cfg <- platform_channels(platform, year)
  out <- list(
    u = as.numeric(data[idx, "chamber_RH"]),
    w = as.matrix(data[idx, cfg$w, drop = FALSE]),
    z = as.matrix(data[idx, cfg$z, drop = FALSE])
  )
  if (references) {
    out$x <- as.numeric(data[idx, cfg$xref])
    out$y <- as.numeric(data[idx, cfg$yref])
  }
  out
}

export_platform_parameters <- function(data, platform, year) {
  cfg <- platform_channels(platform, year)
  d <- complete_platform_data(data, platform, year)
  refs <- c(cfg$xref, cfg$xref, cfg$yref, cfg$yref)
  sensors <- c(cfg$w, cfg$z)
  fit_set <- function(spec) {
    lapply(seq_along(sensors), function(j) {
      fit_calibration(d[[refs[j]]], d[[sensors[j]]], d$chamber_RH, spec)
    })
  }
  fits_u_add <- fit_set("u_add")
  fits_u_gain <- fit_set("u_gain")
  residuals <- matrix(NA_real_, nrow(d), length(sensors))
  for (j in seq_along(sensors)) {
    dp <- calibration_functions(fits_u_add[[j]], d$chamber_RH)
    residuals[, j] <- d[[sensors[j]]] - (dp$delta + dp$phi * d[[refs[j]]])
  }
  list(
    n = nrow(d),
    fits_u_add = fits_u_add,
    fits_u_gain = fits_u_gain,
    raw_sensor_error_cov = cov(residuals),
    raw_sensor_error_residuals = residuals
  )
}

build_epa_parameters <- function(data_list = read_epa_data()) {
  p_a10_21 <- export_platform_parameters(data_list$`2021`, "A10", "2021")
  p_a10_19 <- export_platform_parameters(data_list$`2019`, "A10", "2019")
  p_air_21 <- export_platform_parameters(data_list$`2021`, "AIRPRO", "2021")
  base <- complete_platform_data(data_list$`2021`, "A10", "2021")
  x <- base$Reference_CO; y <- base$Reference_NO2; u <- base$chamber_RH
  mx <- mean(x); sx <- sd(x); my <- mean(y); sy <- sd(y)
  qx <- (x - mx) / sx; qy <- (y - my) / sy
  rh_fit <- lm.fit(cbind(1, qx, qy), u)
  list(
    latent_2021 = list(mean_x = mx, sd_x = sx, mean_y = my, sd_y = sy, rho = cor(x, y)),
    rh_given_latent = list(gamma = rh_fit$coefficients, resid_sd = sd(rh_fit$residuals)),
    epa_2021_latent_rh = as.matrix(base[, c("Reference_CO", "Reference_NO2", "chamber_RH")]),
    platforms = list(
      A10 = list(`2021` = p_a10_21, `2019` = p_a10_19),
      AIRPRO = list(`2021` = p_air_21)
    )
  )
}

raw_channel_coefficients <- function(params, platform, spec = "u_add", campaign = "2021") {
  fits <- params$platforms[[platform]][[campaign]][[paste0("fits_", spec)]]
  convert <- function(fit) {
    b <- fit$beta; u0 <- fit$u0; us <- fit$us
    if (spec == "u_add") return(c(b[1] - b[2] * u0 / us, b[2] / us, b[3]))
    c(b[1] - b[2] * u0 / us, b[2] / us, b[3] - b[4] * u0 / us, b[4] / us)
  }
  t(vapply(fits, convert, numeric(if (spec == "u_add") 3 else 4)))
}

transport_coefficients <- function(params, platform, factor) {
  if (platform != "A10") stop("transport stress is defined for A10")
  c21 <- raw_channel_coefficients(params, platform, "u_add", "2021")
  c19 <- raw_channel_coefficients(params, platform, "u_add", "2019")
  c21 + factor * (c19 - c21)
}
