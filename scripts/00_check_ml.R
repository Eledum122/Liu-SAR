# Check: our concentrated ML must match spatialreg::lagsarlm
suppressMessages({library(spData); library(sf); library(spdep); library(spatialreg)})
source("R/liu_sar_functions.R")

# Read Boston tracts (file is .gpkg or .shp depending on spData version)
f  <- list.files(system.file("shapes", package = "spData"),
                 pattern = "^boston_tracts\\.(gpkg|shp)$", full.names = TRUE)[1]
b  <- st_read(f, quiet = TRUE)

v  <- c("INDUS", "AGE", "DIS", "RAD", "TAX", "CRIM", "LSTAT")  # explanatory variables
lw <- nb2listw(poly2nb(b), style = "W")   # Queen neighbors, row-standardized weights
W  <- listw2mat(lw)                       # weights as a full n x n matrix
y  <- log(b$NOX)                          # response: log of NOX
px <- prep_X(st_drop_geometry(b)[, v])    # center and scale X to unit length

fit <- liu_sar(y, px$X, W)                # our estimators (ML, ridge, Liu)

# Reference ML fit from spatialreg (intercept added automatically;
# "eigen" uses the eigenvalues of W for the log-determinant, as in our code)
ref <- lagsarlm(y ~ px$X, listw = lw, method = "eigen")

# Compare rho, beta, and sigma^2 (coef(ref) = rho, intercept, beta)
cat("rho: ours =", round(fit$rho, 6), " lagsarlm =", round(ref$rho, 6), "\n",
    "max |beta diff| =", signif(max(abs(fit$est[, "ML"] - coef(ref)[-(1:2)])), 3), "\n",
    "sigma2: ours =", round(fit$sigma2, 6), " lagsarlm =", round(ref$s2, 6), "\n")

# Quick sanity run of the remaining outputs (506 tracts, not final)
print(round(fit$d, 3)); cat("k =", round(fit$k, 4), "\n")   # shrinkage parameters
print(round(fit$est, 4))                                    # all estimators side by side
print(lapply(sar_impacts(fit, fit$d["d1"]), round, 4))      # impacts using d1
