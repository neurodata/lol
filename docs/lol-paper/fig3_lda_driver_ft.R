# Parallelize Stuff
#=========================#
library(parallel)

no_cores = detectCores() - 5

cl = makeCluster(no_cores)

# Setup Sims
#==========================#
require(lol)
n=100
niter <- 500  # number of iterations per simulation
rlen <- 30
# the simulations to call themselves
#sims <- list(lol.sims.rtrunk, lol.sims.toep, lol.sims.rtrunk, lol.sims.fat_tails, lol.sims.qdtoep)
#maxr <- c(30, 90, 30, 30, 30)
#ds <- c(100, 100, 100, 1000, 100)
# additional arguments for each simulation scenario
#opt_args <- list(list(), list(), list(K=3), list(rotate=TRUE, priors=c(0.8, 0.2)), list())
#sim_names = c("Trunk-2", "Toeplitz", "Trunk-3", "Fat-Tails (D=1000)", "QDA")
sims <- c(lol.sims.fat_tails)
maxr <- c(30)
ds <- c(1000)
opt_args <- list(list(rotate=TRUE, priors=c(0.8, 0.2)))
sim_names <- c("Fat-Tails (D=1000)")
simulations <- list()
counter <- 1

for (i in 1:length(sims)) {
  for (j in 1:niter) {
    simulations[[counter]] <- list(sim_func=sims[[i]], args=c(list(n, ds[i]), opt_args[[i]]),
                                   rmax=maxr[i], sim=sim_names[i], iter=j)
    counter <- counter + 1
  }
}

# Setup Algorithms
#=========================#
algs <- list(lol.project.pca, lol.project.cpca, lol.project.lrcca, lol.project.lol)
alg_name <- c("PCA", "cPCA", "LR-CCA", "LOL")

clusterExport(cl, "algs"); clusterExport(cl, "alg_name")
clusterExport(cl, "simulations"); clusterExport(cl, "rlen")
results <- parLapply(cl, simulations, function(sim) {
  require(lol)
  sim_dat <- do.call(sim$sim_func, sim$args)
  X <- sim_dat$X; Y <- sim_dat$Y
  results <- data.frame(sim=c(), iter=c(), alg=c(), r=c(), lhat=c())
  for (i in 1:length(algs)) {
    rs <- round(seq(from=1, to=sim$rmax, length.out=rlen))
    for (r in rs) {
      tryCatch({
        xv_res <- lol.xval.eval(X, Y, r, algs[[i]])
        lhat <- xv_res$Lhat
      }, error=function(e) lhat <- NaN)
      results <- rbind(results, data.frame(sim=sim$sim, iter=sim$iter, alg=alg_name[i], r=r, lhat=lhat))
    }
  }
  return(results)
})

# Aggregate and save
#=================================#
require(data.table)
results <- do.call(rbind, results)
results <- data.table(results)
# aggregate over the iterations, retaining the other factors
results.means <- aggregate(lhat ~ sim + alg + r + lhat, data = results, FUN = mean)
results_agg <- list(overall=results, means=results.means)
saveRDS(results_agg, 'lol_fig3_lda_ft.rds')
stopCluster(cl)
