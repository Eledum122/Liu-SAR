# COVID-19 France application (data from Chavez-Chong et al., RRSARMMI repository)
suppressMessages({library(sf); library(spdep)})
source("R/liu_sar_functions.R")

d <- st_read("data/covid/covid_data.shp", quiet = TRUE)
v <- c("LnPop", "C3", "SLivR", "Work", "Inac", "A65pls", "Emer", "FDoc")  # their model
y <- d$LnHosp

# ---- Collinearity ----
R <- cor(st_drop_geometry(d)[, v]); e <- eigen(R)$values
cat("maxVIF =", round(max(diag(solve(R))), 2), " lambda_min =", round(min(e), 4),
    " CN =", round(sqrt(max(e) / min(e)), 2),
    " pairs |r|>0.7 =", sum(abs(R[upper.tri(R)]) > 0.7), "\nVIF:\n")
print(round(diag(solve(R)), 2))

# ---- Two weight matrices ----
nb <- poly2nb(d, queen = TRUE)
# (1) As in Chavez-Chong et al.: inverse-distance queen weights, divided by largest eigenvalue
dl <- nbdists(nb, st_coordinates(st_centroid(st_geometry(d))), longlat = TRUE)
Wcc <- listw2mat(nb2listw(nb, glist = lapply(dl, function(x) 1 / x), style = "B"))
Wcc <- Wcc / max(Re(eigen(Wcc, only.values = TRUE)$values))
# (2) Row-standardized queen weights (as assumed in the paper)
Wrs <- listw2mat(nb2listw(nb, style = "W"))

px <- prep_X(st_drop_geometry(d)[, v])
for (nm in c("Wcc", "Wrs")) {
  fit <- liu_sar(y, px$X, get(nm)); lam <- fit$lambda
  m0  <- est_mse(fit, 1)
  cat("\n==", nm, "== rho =", round(fit$rho, 4), " k =", round(fit$k, 4), "\n")
  print(round(fit$d, 4))
  cat("MSE ratio: Ridge =", round(est_mse(fit, lam / (lam + fit$k)) / m0, 4),
      " Liu(d1) =", round(est_mse(fit, (lam + fit$d[["d1"]]) / (lam + 1)) / m0, 4),
      " Liu(d3) =", round(est_mse(fit, (lam + fit$d[["d3"]]) / (lam + 1)) / m0, 4), "\n")
  print(round(fit$est, 4))
}

# ---- Comparison with published SAR and RRSAR (Chavez-Chong et al., original units) ----
fit <- liu_sar(y, px$X, Wcc)
orig <- fit$est[, c("ML", "Ridge", "d1", "d3")] / px$norm          # back to original units
pub  <- cbind(SAR_pub   = c(-0.1256, 0.5337, 0.6277, 12.8074, 6.7543, -0.0357, -6.6582, -0.1584),
              RRSAR_pub = c(-0.0006, 0.3333, 0.2593,  6.7131, 5.1264, -0.0257, -4.8059, -0.2118))
tab4 <- cbind(orig, pub)
cat("\nOriginal-unit coefficients, Wcc (published: rho_SAR = 0.8489, rho_RRSAR = 0.7282)\n")
print(round(tab4, 4))
write.csv(round(tab4, 6), "output/tables/table4_covid_comparison.csv")

# ---- Diagnostics and summary measures for the paper (Wcc) ----
ols <- lm(y ~ ., data = st_drop_geometry(d)[, v])
lwcc <- mat2listw(Wcc, style = "M")
mi  <- lm.morantest(ols, lwcc, zero.policy = TRUE)
lmt <- if (exists("lm.RStests")) lm.RStests(ols, lwcc, test = "all", zero.policy = TRUE) else
       lm.LMtests(ols, lwcc, test = "all", zero.policy = TRUE)
cat("\nMoran's I =", round(mi$estimate[1], 3), " p =", signif(mi$p.value, 3), "\nLM p-values:",
    paste(names(lmt), signif(sapply(lmt, function(t) t$p.value), 3), sep = "=", collapse = "  "), "\n")
lam <- fit$lambda; m0 <- est_mse(fit, 1)
cat("d:", round(fit$d, 3), "\nMSE ratios: Ridge", round(est_mse(fit, lam/(lam+fit$k))/m0, 3),
    sapply(c("d1","d2","d3"), function(r) round(est_mse(fit, (lam+fit$d[[r]])/(lam+1))/m0, 3)), "\n")
seT <- function(dd) sar_impacts(fit, dd)$Total[, "se"]
cat("a_T =", round(fit$mult$T["a"], 3), " mean SE ratio total (d4T) =",
    round(mean(seT(fit$d[["d4T"]]) / seT(1)), 3), "\n")
print(round(cbind(ML = sar_impacts(fit, 1)$Total, Liu = sar_impacts(fit, fit$d[["d4T"]])$Total) /
            px$norm, 4))   # total impacts, original units
cat("components:", nb_c <- n.comp.nb(nb)$nc, "\n")
