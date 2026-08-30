required <- c(
  "results/simulation/primary_summary.csv",
  "results/simulation/robustness_summary.csv",
  "results/epa/full_sample_results.csv",
  "results/epa/heldout_validation_summary.csv"
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("missing result files: ", paste(missing, collapse = ", "))

prim <- read.csv("results/simulation/primary_summary.csv", stringsAsFactors = FALSE)
stress <- read.csv("results/simulation/robustness_summary.csv", stringsAsFactors = FALSE)
epa <- read.csv("results/epa/full_sample_results.csv", stringsAsFactors = FALSE)
hold <- read.csv("results/epa/heldout_validation_summary.csv", stringsAsFactors = FALSE)
dir.create("figures", showWarnings = FALSE)

cols <- c("#4C78A8", "#F58518", "#54A24B", "#E45756")
pch <- c(16, 15, 17, 18)

a10 <- prim[prim$scenario_id %in% sprintf("P%02d", 0:4), ]
a10 <- a10[match(sprintf("P%02d", 0:4), a10$scenario_id), ]
x <- seq_len(nrow(a10))
labs <- paste0("(", a10$n_collocation, ",", a10$n_target, ")")

png("figures/simulation_interval_performance.png", 2400, 900, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1), oma = c(0, 0, 0, 0), family = "serif")
y <- cbind(a10$jel_coverage, a10$ajel_coverage, a10$bc_jel_coverage, a10$bc_ajel_coverage)
matplot(x, y, type = "b", lty = 1, pch = pch, col = cols, xaxt = "n",
        ylim = c(0.82, 0.985), xlab = "(m, n)", ylab = "Coverage probability",
        main = "(a) 95% confidence-interval coverage")
axis(1, x, labs)
abline(h = .95, lty = 2)
y <- cbind(a10$jel_avg_length, a10$ajel_avg_length, a10$bc_jel_avg_length, a10$bc_ajel_avg_length)
matplot(x, y, type = "b", lty = 1, pch = pch, col = cols, xaxt = "n",
        xlab = "(m, n)", ylab = "Average interval length",
        main = "(b) Average confidence-interval length")
axis(1, x, labs)
legend("topright", c("JEL", "AJEL", "BC-JEL", "BC-AJEL"),
       lty = 1, pch = pch, col = cols, bty = "n")
dev.off()

png("figures/simulation_point_performance.png", 2400, 900, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1), family = "serif")
matplot(x, cbind(a10$proxy_bias, a10$point_bias), type = "b", lty = 1,
        pch = c(16, 17), col = cols[1:2], xaxt = "n", xlab = "(m, n)",
        ylab = "Bias", main = "(a) Bias")
axis(1, x, labs)
abline(h = 0, lty = 2)
legend("bottomright", c("Calibration proxy", "Corrected estimator"),
       lty = 1, pch = c(16, 17), col = cols[1:2], bty = "n")
matplot(x, cbind(a10$proxy_rmse, a10$point_rmse), type = "b", lty = 1,
        pch = c(16, 17), col = cols[1:2], xaxt = "n", xlab = "(m, n)",
        ylab = "RMSE", main = "(b) Root mean squared error")
axis(1, x, labs)
dev.off()

spec_lab <- c(const = "Constant", u_add = "RH additive", u_gain = "RH varying gain")
const <- epa[epa$spec == "const", ]
const <- const[match(c("A10", "AIRPRO"), const$platform), ]
proxy <- hold[hold$spec == "const" & hold$method == "proxy", ]
corr <- hold[hold$spec == "const" & hold$method == "corrected", ]
imp <- hold[hold$spec == "const" & hold$method == "corrected_vs_proxy", ]
proxy <- proxy[match(c("A10", "AIRPRO"), proxy$platform), ]
corr <- corr[match(c("A10", "AIRPRO"), corr$platform), ]
imp <- imp[match(c("A10", "AIRPRO"), imp$platform), ]
rm <- hold[hold$method == "corrected", ]
rmse <- sapply(c("A10", "AIRPRO"), function(platform) {
  g <- rm[rm$platform == platform, ]
  g$rmse[match(c("const", "u_add", "u_gain"), g$spec)]
})
rownames(rmse) <- spec_lab
colnames(rmse) <- c("Kunak A10", "Kunak Air Pro")

png("figures/epa_correlation_validation.png", 2400, 1800, res = 200)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1), family = "serif")
b <- rbind(const$reference_rho, const$epa_study_calibrated_rho, const$corrected_rho)
rownames(b) <- c("Reference", "Study calibrated", "Corrected")
barplot(b, beside = TRUE, names.arg = c("Kunak A10", "Kunak Air Pro"),
        ylim = c(.65, .875), ylab = "Pearson correlation",
        main = "(a) Full-sample correlation diagnostic",
        col = c("grey70", cols[1], cols[2]), density = c(0, 18, 18), angle = c(0, 90, 45))
legend("top", rownames(b), fill = c("grey70", cols[1], cols[2]), bty = "n", horiz = TRUE)

b <- rbind(proxy$mae, corr$mae)
rownames(b) <- c("Calibration proxy", "Corrected estimator")
barplot(b, beside = TRUE, names.arg = c("Kunak A10", "Kunak Air Pro"),
        ylab = "Mean absolute error", main = "(b) Repeated held-out-burn validation",
        col = cols[1:2], density = 18, angle = c(90, 45))
legend("topright", rownames(b), fill = cols[1:2], bty = "n")

barplot(rmse, beside = TRUE, names.arg = c("Kunak A10", "Kunak Air Pro"),
        ylab = "RMSE", main = "(c) Corrected-estimator RMSE by calibration model",
        col = cols[1:3])
legend("topright", rownames(rmse), fill = cols[1:3], bty = "n")

b <- rbind(100 * imp$fraction_improved, 100 * imp$relative_mae_reduction)
rownames(b) <- c("Splits improved (%)", "MAE reduction (%)")
barplot(b, beside = TRUE, names.arg = c("Kunak A10", "Kunak Air Pro"),
        ylim = c(0, 80), ylab = "Percent", main = "(d) Improvement over calibration proxy",
        col = cols[1:2], density = 18, angle = c(90, 45))
legend("topright", rownames(b), fill = cols[1:2], bty = "n")
dev.off()

ids <- c("R01", "R03", "R04", "E02", "D01", "D02", "M01", "M03", "T01", "T02")
s <- stress[match(ids, stress$scenario_id), ]
x <- seq_along(ids)
png("figures/supplement_robustness_performance.png", 1800, 1800, res = 200)
par(mfrow = c(2, 1), mar = c(4.5, 4.5, 2.5, 1), family = "serif")
matplot(x, cbind(s$proxy_bias, s$point_bias), type = "b", lty = 1,
        pch = c(16, 17), col = cols[1:2], xaxt = "n", xlab = "Scenario",
        ylab = "Bias", main = "(a) Point-estimator bias under robustness settings")
axis(1, x, ids)
abline(h = 0, lty = 2)
legend("bottomleft", c("Calibration proxy", "Corrected estimator"),
       lty = 1, pch = c(16, 17), col = cols[1:2], bty = "n")
matplot(x, cbind(s$jel_coverage, s$ajel_coverage, s$bc_ajel_coverage),
        type = "b", lty = 1, pch = c(16, 15, 18), col = cols[c(1, 2, 4)],
        xaxt = "n", ylim = c(.88, .99), xlab = "Scenario",
        ylab = "Coverage probability", main = "(b) 95% confidence-interval coverage")
axis(1, x, ids)
abline(h = .95, lty = 2)
legend("bottomleft", c("JEL", "AJEL", "BC-AJEL"), lty = 1,
       pch = c(16, 15, 18), col = cols[c(1, 2, 4)], bty = "n")
dev.off()
