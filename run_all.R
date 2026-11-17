# Run all analyses in order (set the working directory to this folder, or open Liu-SAR.Rproj)
source("scripts/00_check_ml.R")          # check ML against spatialreg::lagsarlm
source("scripts/boston/01_data.R")       # Boston: 96 NOX zones + diagnostics
source("scripts/boston/02_estimate.R")   # Boston: paper Tables 2-3
source("scripts/boston/03_sensitivity.R")# Boston: paper Table 4
source("scripts/boston/04_figures.R")    # Boston: Figures 1-3
source("scripts/covid/01_estimate.R")    # COVID France: paper Table 5
source("scripts/simulation/02_tables.R")  # Simulation: paper Tables 6-7 (and supplementary Tables 8-9), Figures 4-5 (needs output/sim/results_R5000.rds)
