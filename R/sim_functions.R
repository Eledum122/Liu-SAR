# =============================================================
# Monte Carlo simulation functions (paper Section 6)
# Requires R/liu_sar_functions.R and the spdep package
# =============================================================

# ---- Spatial weights on a sqrt(n) x sqrt(n) grid (row-standardized) ----
make_W <- function(n, type = c("rook", "knn4")) {
  type <- match.arg(type); s <- sqrt(n)
  if (type == "rook") {
    nb <- spdep::cell2nb(s, s, type = "rook")
  } else {                                             # jittered grid avoids ties
    set.seed(1000 + n)
    xy <- as.matrix(expand.grid(1:s, 1:s)) + matrix(runif(2 * n, -0.25, 0.25), n)
    nb <- spdep::knn2nb(spdep::knearneigh(xy, k = 4))
  }
  W <- spdep::listw2mat(spdep::nb2listw(nb, style = "W"))
  list(W = W, om = eigen(W, only.values = TRUE)$values)
}

# ---- Collinear regressors: corr(x_j, x_k) = gamma^2 (Kibria, 2003) ----
make_X <- function(n, p, gamma, seed) {
  set.seed(seed)
  Z <- matrix(rnorm(n * (p + 1)), n)
  sqrt(1 - gamma^2) * Z[, 1:p] + gamma * Z[, p + 1]
}

# ---- One scenario: R replications with fixed X and W ----
# Errors are measured on the original scale of the regressors.
# bdir: true beta along the eigenvector of the largest ("max") or smallest ("min")
#       eigenvalue of R_X (Newhouse and Oman, 1971), scaled so that sum(beta^2) = 1
run_scenario <- function(n, rho, gamma, wtype, sigma, bdir = c("max", "min"),
                         p = 4, R = 1000, alpha = 1, seed = 1) {
  bdir <- match.arg(bdir)
  Wl <- make_W(n, wtype); W <- Wl$W
  Xo <- make_X(n, p, gamma, seed = 100 * n + round(100 * gamma))  # same X for all rho, sigma, W
  px <- prep_X(Xo); X <- px$X
  ev <- eigen(crossprod(X), symmetric = TRUE)$vectors
  bt <- ev[, if (bdir == "max") 1 else p] / px$norm               # true beta (original scale)
  bt <- bt / sqrt(sum(bt^2))
  S  <- solve(diag(n) - rho * W)
  mu <- drop(S %*% (alpha + Xo %*% bt))
  m  <- multipliers(rho, Wl$om, W)
  psi <- sapply(c(D = "D", I = "I", T = "T"), function(k) m[[k]]["a"] * bt)  # true impacts

  dgrid <- seq(0, 1, by = 0.05)                                   # for the oracle fixed d
  est_names <- c("ML", "Ridge", "d1", "d2", "d3", "d4D", "d4I", "d4T")
  set.seed(seed)
  out <- replicate(R, simplify = FALSE, {
    y   <- mu + sigma * drop(S %*% rnorm(n))
    fit <- liu_sar(y, X, W, Wl$om)
    B   <- fit$est / px$norm                                      # original scale
    mh  <- fit$mult
    # squared errors of coefficients
    se_b <- colSums((B - bt)^2)
    # impacts: ML, Ridge, Liu with impact-based d (and d3), original scale
    se_imp <- sapply(c("D", "I", "T"), function(k) {
      a <- mh[[k]]["a"]
      c(ML = sum((a * B[, "ML"] - psi[, k])^2),
        Ridge = sum((a * B[, "Ridge"] - psi[, k])^2),
        d3 = sum((a * B[, "d3"] - psi[, k])^2),
        d4 = sum((a * B[, paste0("d4", k)] - psi[, k])^2))
    })
    # Liu over the d grid (for the oracle)
    BL <- sapply(dgrid, function(d) drop(fit$Fd(d) %*% fit$est[, "ML"]) / px$norm)
    grid_b <- colSums((BL - bt)^2)
    grid_imp <- sapply(c("D", "I", "T"), function(k)
      colSums((mh[[k]]["a"] * BL - psi[, k])^2))
    list(se_b = se_b, se_imp = se_imp, grid_b = grid_b, grid_imp = grid_imp,
         rho = fit$rho, d = fit$d, k = fit$k)
  })

  # ---- Summaries ----
  mse_b   <- rowMeans(sapply(out, `[[`, "se_b"))
  mse_imp <- Reduce(`+`, lapply(out, `[[`, "se_imp")) / R
  g_b     <- rowMeans(sapply(out, `[[`, "grid_b"))
  g_imp   <- Reduce(`+`, lapply(out, `[[`, "grid_imp")) / R
  rh      <- sapply(out, `[[`, "rho")
  res <- c(n = n, rho = rho, gamma = gamma, sigma = sigma, R = R,
    lambda_min = min(eigen(crossprod(X), only.values = TRUE)$values),
    setNames(mse_b, paste0("mse_", names(mse_b))),
    mse_oracle = min(g_b), d_oracle = dgrid[which.min(g_b)],
    setNames(as.vector(mse_imp), paste0("mse_", rep(colnames(mse_imp), each = 4), "_",
                                         rownames(mse_imp))),
    setNames(apply(g_imp, 2, min), paste0("mse_", colnames(g_imp), "_oracle")),
    rho_bias = mean(rh) - rho, rho_rmse = sqrt(mean((rh - rho)^2)),
    setNames(rowMeans(sapply(out, `[[`, "d")), paste0("mean_", names(out[[1]]$d))),
    mean_k = mean(sapply(out, `[[`, "k")))
  data.frame(W = wtype, bdir = bdir, as.list(res))
}
