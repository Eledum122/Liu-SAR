# Boston application - Step 4: estimation and tables (set A, 96 zones)
suppressMessages({library(sf); library(spdep)})
source("R/liu_sar_functions.R")

zones <- readRDS("data/boston_zones.rds")              # from scripts/boston/01_data.R
setA  <- c("INDUS", "AGE", "DIS", "RAD", "TAX", "CRIM", "LSTAT")
W     <- listw2mat(nb2listw(poly2nb(zones), style = "W"))
y     <- log(zones$NOX)
px    <- prep_X(st_drop_geometry(zones)[, setA])
fit   <- liu_sar(y, px$X, W)

# ---- Table 1: coefficients, shrinkage parameters, and estimated MSE ----
lam <- fit$lambda
f_liu <- function(d) (lam + d) / (lam + 1)
mse <- c(ML    = est_mse(fit, 1),
         Ridge = est_mse(fit, lam / (lam + fit$k)),
         sapply(fit$d, function(d) est_mse(fit, f_liu(d))))
tab1 <- rbind(fit$est,
              "shrinkage (k or d)" = c(NA, fit$k, fit$d),
              "estimated MSE"      = mse,
              "MSE / MSE(ML)"      = mse / mse["ML"])
cat("rho =", round(fit$rho, 4), " sigma2 =", signif(fit$sigma2, 4),
    " v_rho =", signif(fit$v_rho, 3), "\n\nTable 1\n")
print(round(tab1, 4))

# ---- Table 2: spatial impacts (ML vs Liu with impact-based d4_kappa) ----
imp_ml <- sar_impacts(fit, 1)
dk     <- c(Direct = "d4D", Indirect = "d4I", Total = "d4T")
tab2 <- do.call(cbind, lapply(names(dk), function(k) {
  liu <- sar_impacts(fit, fit$d[dk[k]])[[k]]
  out <- cbind(imp_ml[[k]], liu)
  colnames(out) <- paste(k, c("ML", "ML.se", "Liu", "Liu.se"))
  out
}))
cat("\nTable 2 (Liu impacts use d4D, d4I, d4T respectively)\n")
print(round(tab2, 4))

write.csv(round(tab1, 6), "output/tables/table1_boston_coef.csv")
write.csv(round(tab2, 6), "output/tables/table2_boston_impacts.csv")
