# Boston application - Step 2: data at 96 NOX zones, Step 3: diagnostics
suppressMessages({library(spData); library(sf); library(spdep); library(dplyr)})
source("R/liu_sar_functions.R")

# ---- Read tracts ----
f <- list.files(system.file("shapes", package = "spData"),
                pattern = "^boston_tracts\\.(gpkg|shp)$", full.names = TRUE)[1]
b <- st_read(f, quiet = TRUE)

# ---- Aggregate 506 tracts to 96 NOX zones ----
# NOX is constant within each zone; covariates are population-weighted means
cand <- c("CRIM", "ZN", "INDUS", "AGE", "DIS", "RAD", "TAX", "PTRATIO", "LSTAT")
zones <- b |>
  group_by(NOX_ID) |>
  summarise(across(all_of(cand), ~ weighted.mean(.x, pmax(POP, 1))),
            NOX = first(NOX), n_tracts = n(), POP = sum(POP),
            do_union = TRUE)
saveRDS(zones, "data/boston_zones.rds")

# ---- Collinearity diagnostics for two candidate sets ----
sets <- list(
  A = c("INDUS", "AGE", "DIS", "RAD", "TAX", "CRIM", "LSTAT"),   # previous set
  B = c("INDUS", "ZN", "AGE", "DIS", "RAD", "TAX")               # land use / urban form
)
collin <- function(X) {
  R <- cor(X); e <- eigen(R)$values
  c(maxVIF = max(diag(solve(R))), lambda_min = min(e),
    CN = sqrt(max(e) / min(e)), n_r_gt_0.7 = sum(abs(R[upper.tri(R)]) > 0.7))
}
Z <- st_drop_geometry(zones)
cat("== Collinearity (96 zones) ==\n")
print(round(t(sapply(sets, function(s) collin(Z[, s]))), 3))
cat("\nVIF, set B:\n"); print(round(diag(solve(cor(Z[, sets$B]))), 2))

# ---- Spatial weights and spatial dependence tests ----
nb <- poly2nb(zones); lw <- nb2listw(nb, style = "W")
cat("\nNeighbors per zone:", summary(card(nb)), "\n")
y <- log(zones$NOX)
for (s in names(sets)) {
  ols <- lm(y ~ ., data = Z[, sets[[s]]])
  mi  <- lm.morantest(ols, lw)
  lm_t <- if (exists("lm.RStests")) lm.RStests(ols, lw, test = "all")
          else lm.LMtests(ols, lw, test = "all")
  pv  <- sapply(lm_t, function(t) t$p.value)
  cat("\nSet", s, ": Moran's I (OLS residuals) =", round(mi$estimate[1], 3),
      " p =", signif(mi$p.value, 3), "\n  LM p-values:",
      paste(names(pv), signif(pv, 3), sep = "=", collapse = "  "), "\n")
}
