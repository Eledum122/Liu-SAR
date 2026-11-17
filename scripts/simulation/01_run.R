# Monte Carlo simulation (Section 6): runs all scenarios in parallel
# Each scenario is saved separately in output/sim/, so the run can be resumed.
suppressMessages(library(parallel))
dir.create("output/sim", showWarnings = FALSE, recursive = TRUE)
R_REP <- as.integer(Sys.getenv("R_REP", "1000"))   # set in run_sim.R

grid <- expand.grid(n = c(49, 100, 225, 400), rho = c(0.2, 0.5, 0.8),
                    gamma = c(0.8, 0.9, 0.95, 0.99), wtype = c("rook", "knn4"),
                    sigma = c(0.5, 1, 2), bdir = c("max", "min"), stringsAsFactors = FALSE)
grid$file <- with(grid, sprintf("output/sim/n%d_rho%.1f_g%.2f_%s_s%.1f_b%s_R%d.rds",
                                n, rho, gamma, wtype, sigma, bdir, R_REP))
todo <- which(!file.exists(grid$file))
cat(nrow(grid), "scenarios,", length(todo), "to run\n")
t0 <- Sys.time()

cl <- makeCluster(as.integer(Sys.getenv("N_CORES", max(1, detectCores() - 1))))    # PSOCK cluster: works on Windows and Linux
invisible(clusterEvalQ(cl, {suppressMessages(library(spdep))
                  source("R/liu_sar_functions.R"); source("R/sim_functions.R")}))
clusterExport(cl, c("grid", "R_REP"))
invisible(parLapplyLB(cl, todo, function(i) {
  g <- grid[i, ]
  res <- run_scenario(g$n, g$rho, g$gamma, g$wtype, g$sigma, g$bdir, R = R_REP, seed = i)
  saveRDS(res, g$file)
}))
stopCluster(cl)

res <- do.call(rbind, lapply(grid$file, readRDS))
saveRDS(res, sprintf("output/sim/results_R%d.rds", R_REP))
el <- difftime(Sys.time(), t0, units = "hours")
cat(sprintf("Done: %d scenarios. Elapsed: %.2f hours\n", nrow(res), as.numeric(el)))
writeLines(sprintf("R = %d, %d scenarios run, %.2f hours, finished %s",
                   R_REP, length(todo), as.numeric(el), format(Sys.time())),
           sprintf("output/sim/runtime_R%d.txt", R_REP))
