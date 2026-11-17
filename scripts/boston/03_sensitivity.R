# Boston application - sensitivity: variable set B and the 506-tract level
suppressMessages({library(spData); library(sf); library(spdep)})
source("R/liu_sar_functions.R")

zones <- readRDS("data/boston_zones.rds")
f <- list.files(system.file("shapes", package = "spData"),
                pattern = "^boston_tracts\\.(gpkg|shp)$", full.names = TRUE)[1]
tracts <- st_read(f, quiet = TRUE)

setA <- c("INDUS", "AGE", "DIS", "RAD", "TAX", "CRIM", "LSTAT")
setB <- c("INDUS", "ZN", "AGE", "DIS", "RAD", "TAX")

# One summary row per configuration
summ <- function(data, vars) {
  W   <- listw2mat(nb2listw(poly2nb(data), style = "W"))
  px  <- prep_X(st_drop_geometry(data)[, vars])
  fit <- liu_sar(log(data$NOX), px$X, W)
  lam <- fit$lambda; fL <- function(d) (lam + d) / (lam + 1)
  m0  <- est_mse(fit, 1)
  seT <- function(d) sar_impacts(fit, d)$Total[, "se"]       # total-impact SEs
  c(n = nrow(data), rho = fit$rho, maxVIF = max(diag(solve(cor(px$X)))),
    lambda_min = min(lam), k = fit$k, d1 = fit$d[["d1"]], d3 = fit$d[["d3"]],
    MSEratio_Ridge = est_mse(fit, lam / (lam + fit$k)) / m0,
    MSEratio_d1 = est_mse(fit, fL(fit$d[["d1"]])) / m0,
    MSEratio_d3 = est_mse(fit, fL(fit$d[["d3"]])) / m0,
    SEratio_Total_d4T = mean(seT(fit$d[["d4T"]]) / seT(1)))
}

tab3 <- rbind("A, 96 zones (main)" = summ(zones, setA),
              "B, 96 zones"        = summ(zones, setB),
              "A, 506 tracts"      = summ(tracts, setA))
print(round(t(tab3), 4))
write.csv(round(tab3, 6), "output/tables/table3_boston_sensitivity.csv")
