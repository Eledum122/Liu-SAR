# Estimate the total running time of the simulation on this computer
suppressMessages(library(spdep))
source("R/liu_sar_functions.R"); source("R/sim_functions.R")
R_REP   <- as.integer(Sys.getenv("R_REP", "1000"))                 # set in run_sim.R
n_cores <- as.integer(Sys.getenv("N_CORES", max(1, parallel::detectCores() - 1)))
reps    <- 20                                       # replications used for timing

sec <- sapply(c(49, 100, 225, 400), function(n)
  system.time(run_scenario(n, 0.5, 0.95, "rook", 1, "max", R = reps))["elapsed"] / reps)
names(sec) <- paste0("n=", c(49, 100, 225, 400))
total_h <- 144 * R_REP * sum(sec) / 3600 / n_cores  # 144 scenarios per sample size

cat("Seconds per replication:\n"); print(round(sec, 3))
cat(sprintf("\nWorkers: %d   Estimated total time for R = %d: %.1f hours\n",
            n_cores, R_REP, total_h))
