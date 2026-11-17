# Liu-SAR: A Liu estimator for spatial autoregressive models under multicollinearity

R code that reproduces all results of the paper

> H. Eledum, *A Liu estimator for spatial autoregressive models under multicollinearity: MSE properties, spatial impacts, and application to air pollution data*, submitted to AIMS Mathematics.

## Requirements

R (version 4.3 or later; the paper used R 4.6.1) with the packages
`spData`, `sf`, `spdep`, `spatialreg`, and `dplyr`. The simulation also uses the base package `parallel`.

```r
install.packages(c("spData", "sf", "spdep", "spatialreg", "dplyr"))
```

## Reproducing the results

Set the working directory to this folder (or open `Liu-SAR.Rproj`) and run

```r
source("run_all.R")
```

This reproduces the empirical applications and the simulation tables and figures in a few minutes. It uses the stored simulation results in `output/sim/results_R5000.rds`.

To rerun the full Monte Carlo study (576 scenarios, R = 5000 replications; about 22 hours with seven parallel processes), run

```r
source("scripts/simulation/run_sim.R")
```

The number of replications and of parallel workers are set at the top of `run_sim.R`.

## Contents

| Path | Contents | Paper |
|---|---|---|
| `R/liu_sar_functions.R` | ML–SAR, ridge–SAR, and Liu–SAR estimation; rules for selecting d; impact multipliers and standard errors | Sections 2–4 |
| `R/sim_functions.R` | Weights matrices, regressors, and one simulation scenario | Section 6 |
| `scripts/00_check_ml.R` | Check of the ML estimates against `spatialreg::lagsarlm` | – |
| `scripts/boston/` | Boston NOX application: data at the level of the 96 NOX zones, estimation, sensitivity, figures | Section 5.1–5.3, Tables 2–4, Figures 1–3 |
| `scripts/covid/01_estimate.R` | COVID-19 France application | Section 5.4, Table 5 |
| `scripts/simulation/` | Monte Carlo study: running time, simulation, tables and figures | Section 6, Tables 6–7, Figures 4–5 |
| `data/` | `boston_zones.rds` (created by `scripts/boston/01_data.R`); `covid/` COVID-19 France data | – |
| `output/` | Tables (CSV), figures (PDF), and stored simulation results | – |

Table numbers in the file names of `output/tables/` follow the scripts; the column "Paper" above gives the corresponding tables in the paper.

## Data sources

- Boston data: Harrison and Rubinfeld (1978), corrected by Gilley and Pace (1996), distributed with the R package `spData`.
- COVID-19 France data: Chavez-Chong, Hardouin and Fermin (2025), *AStA Advances in Statistical Analysis* 109, 25–52, from the repository <https://github.com/c0ra/RRSARMMI> (MIT license).

## License

The code is released under the MIT License (see `LICENSE`).
