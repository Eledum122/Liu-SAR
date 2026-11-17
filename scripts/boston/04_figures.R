# Boston application - Step 5: figures (PDF for the paper, PNG for preview)
suppressMessages({library(spData); library(sf); library(spdep)})
source("R/liu_sar_functions.R")
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# ---- Data and fit (same as 02_estimate.R) ----
zones <- readRDS("data/boston_zones.rds")
f <- list.files(system.file("shapes", package = "spData"),
                pattern = "^boston_tracts\\.(gpkg|shp)$", full.names = TRUE)[1]
tracts <- st_read(f, quiet = TRUE)
setA <- c("INDUS", "AGE", "DIS", "RAD", "TAX", "CRIM", "LSTAT")
W   <- listw2mat(nb2listw(poly2nb(zones), style = "W"))
px  <- prep_X(st_drop_geometry(zones)[, setA])
fit <- liu_sar(log(zones$NOX), px$X, W)

# ---- Common style ----
cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6")
blues <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
ink <- "#3d3d3a"; grid_col <- "#e6e5e1"
save_fig <- function(name, w, h, draw) {           # same figure as PDF and PNG
  pdf(file.path("output/figures", paste0(name, ".pdf")), width = w, height = h); draw(); dev.off()
  png(file.path("output/figures", paste0(name, ".png")), width = w, height = h,
      units = "in", res = 200); draw(); dev.off()
}
base_par <- function() par(family = "sans", col.axis = ink, col.lab = ink, fg = ink,
                           las = 1, tcl = -0.25, mgp = c(2.2, 0.5, 0), cex.axis = 0.85)

# ---- Figure 1: log(NOX) for 506 tracts and 96 NOX zones ----
brks <- quantile(log(tracts$NOX), seq(0, 1, length.out = 8))
save_fig("fig1_boston_map", 9, 4.6, function() {
  par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5), oma = c(3, 0, 0, 0))
  plot(st_geometry(tracts), col = blues[cut(log(tracts$NOX), brks, include.lowest = TRUE)],
       border = "white", lwd = 0.2)
  title("(a) 506 census tracts", col.main = ink, font.main = 1, cex.main = 1)
  plot(st_geometry(zones), col = blues[cut(log(zones$NOX), brks, include.lowest = TRUE)],
       border = "white", lwd = 0.6)
  title("(b) 96 NOX zones", col.main = ink, font.main = 1, cex.main = 1)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("bottom", legend = sprintf("%.2f to %.2f", head(brks, -1), brks[-1]),
         fill = blues, border = NA, horiz = TRUE, bty = "n", cex = 0.7, text.col = ink,
         title = "log(NOX)", title.col = ink)
})

# ---- Figure 2: Liu trace and estimated MSE as functions of d ----
dg   <- seq(0, 1, by = 0.01)
lam  <- fit$lambda
B    <- sapply(dg, function(d) drop(fit$Fd(d) %*% fit$est[, "ML"]))     # p x length(dg)
mse  <- sapply(dg, function(d) est_mse(fit, (lam + d) / (lam + 1)))
mR   <- est_mse(fit, lam / (lam + fit$k))
dsel <- fit$d[c("d3", "d4I")]                       # d1 = 0.727 is almost equal to d3
mark_d <- function(ylim) {
  abline(v = dsel, lty = 2, col = "#9a9893", lwd = 1)
  text(dsel, ylim[2], c(expression(hat(d)[3] %~~% hat(d)[1]), expression(hat(d)["4,I"])),
       pos = c(4, 2), cex = 0.8, col = "#6b6a65", xpd = NA)
}
save_fig("fig2_boston_liu_trace_mse", 9, 4.2, function() {
  par(mfrow = c(1, 2), mar = c(3.6, 3.8, 1.8, 4.5)); base_par()
  # (a) trace
  yl <- range(B) * c(1.05, 1.08)
  plot(NA, xlim = c(0, 1), ylim = yl, xlab = "d", ylab = "Coefficient (unit-length scale)",
       panel.first = grid(nx = NA, ny = NULL, col = grid_col, lty = 1), bty = "n")
  abline(h = 0, col = "#c9c8c3")
  mark_d(yl)
  for (j in seq_along(setA)) lines(dg, B[j, ], col = cols[j], lwd = 2)
  yr <- B[, length(dg)]; yo <- order(yr); ylab_pos <- yr
  gap <- diff(yl) * 0.045                                        # avoid label overlap
  for (i in 2:length(yo)) ylab_pos[yo[i]] <- max(yr[yo[i]], ylab_pos[yo[i - 1]] + gap)
  text(1.02, ylab_pos, setA, pos = 4, cex = 0.7, col = ink, xpd = NA)
  title("(a) Liu trace", adj = 0, font.main = 1, cex.main = 1, col.main = ink)
  # (b) estimated MSE
  yl <- c(min(mse, mR) * 0.97, max(mse) * 1.03)
  plot(dg, mse, type = "n", ylim = yl, xlab = "d", ylab = "Estimated MSE",
       panel.first = grid(nx = NA, ny = NULL, col = grid_col, lty = 1), bty = "n")
  mark_d(yl)
  abline(h = mR, lty = 3, col = cols[2], lwd = 1.5)
  text(0, mR, "Ridge-SAR", adj = c(0, -0.6), cex = 0.75, col = ink)
  lines(dg, mse, col = cols[1], lwd = 2)
  points(1, mse[length(dg)], pch = 21, bg = cols[1], col = "white", cex = 1.4)
  text(1, mse[length(dg)], "ML-SAR", pos = 2, cex = 0.75, col = ink)
  title("(b) Estimated MSE of Liu-SAR", adj = 0, font.main = 1, cex.main = 1, col.main = ink)
})

# ---- Figure 3: impacts with 95% intervals, ML-SAR vs Liu-SAR ----
imp_ml <- sar_impacts(fit, 1)
dk <- c(Direct = "d4D", Indirect = "d4I", Total = "d4T")
save_fig("fig3_boston_impacts", 9, 3.8, function() {
  par(mfrow = c(1, 3), mar = c(3.4, 0.6, 1.8, 0.6), oma = c(1.6, 4.2, 0, 0)); base_par()
  y <- rev(seq_along(setA))
  for (k in names(dk)) {
    liu <- sar_impacts(fit, fit$d[dk[k]])[[k]]; ml <- imp_ml[[k]]
    lo <- c(ml[, 1] - 1.96 * ml[, 2], liu[, 1] - 1.96 * liu[, 2])
    hi <- c(ml[, 1] + 1.96 * ml[, 2], liu[, 1] + 1.96 * liu[, 2])
    plot(NA, xlim = range(lo, hi), ylim = c(0.5, length(y) + 0.5), yaxt = "n",
         xlab = "Impact", ylab = "", bty = "n",
         panel.first = grid(nx = NULL, ny = NA, col = grid_col, lty = 1))
    abline(v = 0, col = "#9a9893")
    if (k == "Direct") axis(2, at = y, labels = setA, tick = FALSE, cex.axis = 0.85)
    for (s in 1:2) {
      est <- if (s == 1) ml else liu; yy <- y + c(0.15, -0.15)[s]
      segments(est[, 1] - 1.96 * est[, 2], yy, est[, 1] + 1.96 * est[, 2], yy,
               col = cols[s], lwd = 2)
      points(est[, 1], yy, pch = 21, bg = cols[s], col = "white", cex = 1.3)
    }
    title(k, adj = 0, font.main = 1, cex.main = 1, col.main = ink)
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
  legend("bottom", c("ML-SAR", "Liu-SAR (impact-based d)"), col = cols[1:2], lwd = 2,
         pch = 21, pt.bg = cols[1:2], horiz = TRUE, bty = "n", cex = 0.85, text.col = ink)
})
cat("Figures saved in output/figures\n")
