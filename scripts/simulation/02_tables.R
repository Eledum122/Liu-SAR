# Simulation results: tables (CSV) and figures (PDF/PNG) for Section 6
R_REP <- as.integer(Sys.getenv("R_REP", "5000"))
r <- readRDS(sprintf("output/sim/results_R%d.rds", R_REP))
dir.create("output/tables", showWarnings = FALSE); dir.create("output/figures", showWarnings = FALSE)

RE  <- function(e, k = "") r[[paste0("mse_", k, "ML")]] / r[[paste0("mse_", k, e)]]   # MSE(ML)/MSE(e)
est <- c("Ridge", "d1", "d2", "d3", "oracle")
re  <- sapply(est, RE)
by_factor <- function(f) do.call(rbind, lapply(c("max", "min"), function(b) {
  i <- r$bdir == b
  a <- aggregate(re[i, ], list(level = r[[f]][i]), mean)
  cbind(beta = b, factor = f, a)
}))

# ---- Table 6: coefficient RE by gamma, n, sigma, rho (beta direction panels) ----
tab6 <- do.call(rbind, lapply(c("gamma", "n", "sigma", "rho"), by_factor))
summ <- do.call(rbind, lapply(c("max", "min"), function(b) {
  i <- r$bdir == b
  rbind(cbind(beta = b, factor = "overall", level = "mean", as.data.frame(t(colMeans(re[i, ])))),
        cbind(beta = b, factor = "overall", level = "share RE>1", as.data.frame(t(colMeans(re[i, ] > 1)))),
        cbind(beta = b, factor = "overall", level = "min", as.data.frame(t(apply(re[i, ], 2, min)))))
}))
tab6 <- rbind(tab6, summ)
write.csv(cbind(tab6[, 1:3], round(tab6[, -(1:3)], 3)), "output/tables/table6_sim_coef.csv", row.names = FALSE)

# ---- Table 7: impacts RE (Ridge, d3, impact-based d4, oracle) ----
tab7 <- do.call(rbind, lapply(c("max", "min"), function(b) {
  i <- r$bdir == b
  t(sapply(c(Direct = "D", Indirect = "I", Total = "T"), function(k) {
    p <- paste0(k, "_")
    c(Ridge = mean(RE("Ridge", p)[i]), d3 = mean(RE("d3", p)[i]), d4 = mean(RE("d4", p)[i]),
      oracle = mean(RE("oracle", p)[i]),
      share_d4_better_d3 = mean(r[[paste0("mse_", p, "d4")]][i] < r[[paste0("mse_", p, "d3")]][i]),
      min_d4 = min(RE("d4", p)[i]))
  })) |> as.data.frame() |> cbind(beta = b, impact = c("Direct", "Indirect", "Total"))
}))
write.csv(cbind(tab7[, c("beta", "impact")], round(tab7[, 1:6], 3)),
          "output/tables/table7_sim_impacts.csv", row.names = FALSE)

# ---- Table 8: bias and RMSE of rho-hat; Table 9: selected shrinkage parameters ----
tab8 <- aggregate(r[, c("rho_bias", "rho_rmse")], list(n = r$n, rho = r$rho), mean)
write.csv(round(tab8, 4), "output/tables/table8_sim_rho.csv", row.names = FALSE)
tab9 <- aggregate(r[, c("mean_d1", "mean_d2", "mean_d3", "d_oracle", "mean_k")],
                  list(beta = r$bdir, gamma = r$gamma), mean)
write.csv(cbind(tab9[, 1:2], round(tab9[, -(1:2)], 3)), "output/tables/table9_sim_d.csv", row.names = FALSE)

# ---- Figures ----
cols <- c(Ridge = "#2a78d6", d2 = "#eb6834", d3 = "#1baf7a", oracle = "#6b6a65")
labs <- c(Ridge = "Ridge-SAR", d2 = "Liu-SAR (d2)", d3 = "Liu-SAR (d3)", oracle = "Liu-SAR (oracle d)")
ink <- "#3d3d3a"; grid_col <- "#e6e5e1"
save_fig <- function(name, w, h, draw) {
  pdf(file.path("output/figures", paste0(name, ".pdf")), width = w, height = h); draw(); dev.off()
  png(file.path("output/figures", paste0(name, ".png")), width = w, height = h, units = "in", res = 200)
  draw(); dev.off()
}
base_par <- function() par(col.axis = ink, col.lab = ink, fg = ink, las = 1, tcl = -0.25,
                           mgp = c(2.2, 0.5, 0), cex.axis = 0.85)
g <- sort(unique(r$gamma))

# Figure 4: RE against gamma, by beta direction (log scale, oracle as reference)
save_fig("fig4_sim_re_gamma", 9, 4, function() {
  par(mfrow = c(1, 2), mar = c(3.6, 3.8, 1.8, 0.8), oma = c(2, 0, 0, 0)); base_par()
  for (b in c("max", "min")) {
    i <- r$bdir == b
    m <- sapply(names(cols), function(e) tapply(re[i, e], r$gamma[i], mean))
    plot(NA, xlim = range(g), ylim = range(c(0.95, m)), log = "y", xaxt = "n",
         xlab = expression(gamma), ylab = "Relative efficiency (log scale)", bty = "n",
         panel.first = grid(nx = NA, ny = NULL, col = grid_col, lty = 1))
    axis(1, at = g)
    abline(h = 1, col = "#9a9893", lty = 2)
    for (e in names(cols)) {
      lines(g, m[, e], col = cols[e], lwd = 2, lty = if (e == "oracle") 3 else 1)
      points(g, m[, e], pch = 21, bg = cols[e], col = "white", cex = 1.2)
    }
    title(if (b == "max") "(a) Favorable orientation" else "(b) Unfavorable orientation",
          adj = 0, font.main = 1, cex.main = 1, col.main = ink)
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
  legend("bottom", labs, col = cols, lwd = 2, lty = c(1, 1, 1, 3), pch = 21, pt.bg = cols,
         horiz = TRUE, bty = "n", cex = 0.8, text.col = ink)
})

# Figure 5: average selected d versus the oracle d, by gamma and beta direction
# d1 is omitted because it is almost identical to d3
dc <- c(mean_d2 = "#eb6834", mean_d3 = "#1baf7a", d_oracle = "#6b6a65")
dl <- c("d2", "d3", "oracle d")
save_fig("fig5_sim_d", 9, 4, function() {
  par(mfrow = c(1, 2), mar = c(3.6, 3.8, 1.8, 0.8), oma = c(2, 0, 0, 0)); base_par()
  for (b in c("max", "min")) {
    t9 <- tab9[tab9$beta == b, ]
    plot(NA, xlim = range(g), ylim = c(0, 1), xaxt = "n", xlab = expression(gamma),
         ylab = "Average d", bty = "n", panel.first = grid(nx = NA, ny = NULL, col = grid_col, lty = 1))
    axis(1, at = g)
    for (v in names(dc)) {
      lines(t9$gamma, t9[[v]], col = dc[v], lwd = 2, lty = if (v == "d_oracle") 3 else 1)
      points(t9$gamma, t9[[v]], pch = 21, bg = dc[v], col = "white", cex = 1.2)
    }
    title(if (b == "max") "(a) Favorable orientation" else "(b) Unfavorable orientation",
          adj = 0, font.main = 1, cex.main = 1, col.main = ink)
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
  legend("bottom", dl, col = dc, lwd = 2, lty = c(1, 1, 3), pch = 21, pt.bg = dc,
         horiz = TRUE, bty = "n", cex = 0.85, text.col = ink)
})
cat("Simulation tables and figures saved\n")
