# Parallelize Stuff
#=========================#
library(parallel)

no_cores = detectCores() - 1

cl = makeCluster(no_cores)

# Setup Sims
#==========================#
require(fselect)
n=100; p=100
rs <- seq(1, 30, by=1)
niter <- 100  # number of iterations per simulation
# the simulations to call themselves
sims <- c(fs.sims.trunk, fs.sims.trunk, fs.sims.toep, fs.sims.fat_tails, fs.sims.qdtoep)
# additional arguments for each simulation scenario
opt_args <- list(list(), list(C=3), list(), list(), list())
