# Run the Monte Carlo simulation (start as an RStudio Background Job or with source())
R_REP          <- 5000                                 # number of replications
N_CORES        <- max(1, parallel::detectCores() - 1)  # parallel workers
ESTIMATE_FIRST <- TRUE                                 # print the estimated running time first
RUN            <- TRUE                                 # FALSE = only estimate the time

Sys.setenv(R_REP = R_REP, N_CORES = N_CORES)
if (ESTIMATE_FIRST) source("scripts/simulation/00_time_estimate.R")
if (RUN)            source("scripts/simulation/01_run.R")
