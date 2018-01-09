# Parallelize Stuff
#=========================#
library(parallel)

no_cores = detectCores() - 1

cl = makeCluster(no_cores)

# Setup Sims
#==========================#
require(fselect)
n=100
niter <- 100  # number of iterations per simulation
rlen <- 30
# the simulations to call themselves
sims <- list(fs.sims.rtrunk, fs.sims.toep, fs.sims.rtrunk, fs.sims.fat_tails, fs.sims.qdtoep)
maxr <- c(30, 90, 30, 30, 30)
ds <- c(100, 100, 100, 1000, 100)
# additional arguments for each simulation scenario
opt_args <- list(list(), list(), list(C=3), list(), list())
sim_names = c("Trunk-2", "Toeplitz", "Trunk-3", "Fat-Tails (D=1000)", "QDA")

simulations <- list()
counter <- 1

for (i in 1:length(sims)) {
  for (j in 1:niter) {
    sim_dat <- do.call(sims[[i]], c(list(n, ds[i]), opt_args[[i]]))
    simulations[[counter]] <- list(X=sim_dat$X, Y=sim_dat$Y, rmax=maxr[i], sim=sim_names[i], iter=j)
    counter <- counter + 1
  }
}

# Setup Algorithms
#=========================#
algs <- list(fs.project.pca, fs.project.cpca, fs.project.lrcca, fs.project.lol)
alg_name <- c("PCA", "cPCA", "LR-CCA", "LOL")

clusterExport(cl, "algs"); clusterExport(cl, "alg_name"); clusterExport(cl, "sim")
clusterExport(cl, "rs"); clusterExport(cl, "maxr")
results <- lapply(simulations, function(sim) {
  require(fselect)
  results <- data.frame(sim=c(), iter=c(), alg=c(), r=c(), lhat=c())
  for (i in 1:length(algs)) {
    rs <- round(seq(from=1, to=sim$rmax, length.out=rlen))
    for (r in rs) {
      tryCatch({
        xv_res <- fs.xval.eval(sim$X, sim$Y, r, algs[[i]])
        lhat <- xv_res$Lhat
      }, error=function(e) lhat <- NaN)
      results <- rbind(results, data.frame(sim=sim$sim, iter=sim$iter, alg=alg_name[i], r=r, lhat=lhat))
    }
  }
  return(results)
})

saveRDS(results, 'lol_fig3.rds')
stopCluster(cl)
