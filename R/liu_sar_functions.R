# =============================================================
# Two-stage Liu-SAR estimation (paper Sections 2-4)
# Notation follows the paper: C = X'X = R_X (unit-length scaling)
# =============================================================

# ---- 1. Center and scale X to unit length (Section 2.1) ----
prep_X <- function(X) {
  X  <- as.matrix(X)
  mu <- colMeans(X)
  Xc <- sweep(X, 2, mu)
  nrm <- sqrt(colSums(Xc^2))
  list(X = sweep(Xc, 2, nrm, "/"), mean = mu, norm = nrm)
}

# ---- 2. Concentrated ML for the SAR model (eqs. 3-6) ----
sar_ml <- function(y, X, W, om = NULL) {        # om: eigenvalues of W (optional, for speed)
  n  <- length(y)
  if (is.null(om)) om <- eigen(W, only.values = TRUE)$values
  Ci <- solve(crossprod(X))
  Wy <- drop(W %*% y)
  yc <- y - mean(y); Wyc <- Wy - mean(Wy)          # M y and M W y
  b0 <- drop(Ci %*% crossprod(X, yc))
  bW <- drop(Ci %*% crossprod(X, Wyc))
  e0 <- yc - drop(X %*% b0); eW <- Wyc - drop(X %*% bW)
  lc <- function(r) -n / 2 * log(sum((e0 - r * eW)^2) / n) +
                    sum(log(Mod(1 - r * om)))
  rho <- optimize(lc, c(-0.99, 0.99), maximum = TRUE, tol = 1e-10)$maximum
  list(rho = rho, beta = b0 - rho * bW, bW = bW,
       sigma2 = sum((e0 - rho * eW)^2) / n,
       alpha = mean(y) - rho * mean(Wy), omega = om, Ci = Ci)
}

# ---- 3. Impact multipliers a_k(rho) and derivatives (Section 4.1) ----
# a_T is computed from S = A(rho)^{-1}; it equals 1/(1-rho) for row-standardized W
multipliers <- function(rho, om, W) {
  aD  <- mean(Re(1 / (1 - rho * om)))
  aD1 <- mean(Re(om / (1 - rho * om)^2))
  S1  <- solve(diag(nrow(W)) - rho * W, rep(1, nrow(W)))       # S 1
  aT  <- mean(S1)
  aT1 <- mean(solve(diag(nrow(W)) - rho * W, W %*% S1))        # 1'SWS1 / n
  list(D = c(a = aD, a1 = aD1), T = c(a = aT, a1 = aT1),
       I = c(a = aT - aD, a1 = aT1 - aD1))
}

# ---- 4. Common d-rule: ratio in eq. (dimp) with canonical variances q ----
d_rule <- function(theta2, q, lam) {
  d0 <- sum((theta2 - lam * q) / (lam + 1)^2) / sum((theta2 + q) / (lam + 1)^2)
  max(0, d0)
}

# ---- 5. Main function: ML, ridge-SAR, and Liu-SAR with all d rules ----
liu_sar <- function(y, X, W, om = NULL) {
  n <- length(y); p <- ncol(X)
  ml <- sar_ml(y, X, W, om)
  C  <- crossprod(X); eg <- eigen(C, symmetric = TRUE)
  Q  <- eg$vectors; lam <- eg$values
  b  <- ml$beta; s2 <- ml$sigma2; rho <- ml$rho

  # v_rho from the efficient information J_rho (eq. vrho)
  H   <- W %*% solve(diag(n) - rho * W)
  mu  <- ml$alpha + drop(X %*% b)                     # E(A y) = alpha 1 + X beta
  HXb <- drop(H %*% mu); HXb <- HXb - mean(HXb)
  MxH <- HXb - drop(X %*% (ml$Ci %*% crossprod(X, HXb)))
  J   <- sum(HXb * MxH) / s2 + sum(H * H) + sum(diag(H %*% H)) -
         2 / n * sum(diag(H))^2
  v   <- 1 / J
  bWs <- drop(ml$Ci %*% crossprod(X, H %*% mu))        # beta_W*

  # canonical quantities
  th <- drop(crossprod(Q, b)); eta <- drop(crossprod(Q, bWs))
  q0 <- s2 / lam; m <- multipliers(rho, ml$omega, W)

  d <- c(d1 = d_rule(th^2, q0, lam),
         d2 = d_rule(pmax(th^2 - q0, 0), q0, lam),
         d3 = d_rule(th^2, q0 + v * eta^2, lam),
         sapply(c(D = "D", I = "I", T = "T"), function(k)
           d_rule(th^2, q0 + v * (eta - m[[k]]["a1"] / m[[k]]["a"] * th)^2, lam)))
  names(d)[4:6] <- paste0("d4", c("D", "I", "T"))

  Fd  <- function(dd) solve(C + diag(p), C + dd * diag(p))
  k   <- p * s2 / sum(b^2)                                  # HKB rule
  est <- cbind(ML = b, Ridge = drop(solve(C + k * diag(p), C %*% b)),
               sapply(d, function(dd) drop(Fd(dd) %*% b)))
  rownames(est) <- colnames(X)

  list(est = est, d = d, k = k, rho = rho, sigma2 = s2, v_rho = v,
       bWs = bWs, mult = m, C = C, Ci = ml$Ci, Fd = Fd, lambda = lam,
       ml = ml)                                             # full sar_ml output
}

# ---- 6. Spatial impacts with delta-method SEs (eq. impvar) ----
sar_impacts <- function(fit, d) {
  Fm <- fit$Fd(d); b <- drop(Fm %*% fit$est[, "ML"])
  out <- lapply(c("D", "I", "T"), function(k) {
    a <- fit$mult[[k]]["a"]; a1 <- fit$mult[[k]]["a1"]
    g <- drop(Fm %*% (a1 * fit$est[, "ML"] - a * fit$bWs))
    V <- a^2 * fit$sigma2 * Fm %*% fit$Ci %*% Fm + fit$v_rho * tcrossprod(g)
    cbind(impact = a * b, se = sqrt(diag(V)))
  })
  names(out) <- c("Direct", "Indirect", "Total")
  out
}

# ---- 7. Estimated scalar MSE (plug-in; includes v_rho term, Section 3.5) ----
# f = canonical shrinkage factors: ML f = 1, Liu f = (lam+d)/(lam+1),
# ridge f = lam/(lam+k).  MSE = sum f^2 tau^2 + (1-f)^2 theta^2
est_mse <- function(fit, f) {
  Q   <- eigen(fit$C, symmetric = TRUE)$vectors
  th  <- drop(crossprod(Q, fit$est[, "ML"]))
  eta <- drop(crossprod(Q, fit$bWs))
  tau2 <- fit$sigma2 / fit$lambda + fit$v_rho * eta^2
  sum(f^2 * tau2 + (1 - f)^2 * th^2)
}
