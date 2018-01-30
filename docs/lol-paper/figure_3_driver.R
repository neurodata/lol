# Parallelize Stuff
#=========================#
library(parallel)

no_cores = detectCores() - 5

cl = makeCluster(no_cores)

# Setup Sims
#==========================#
require(lol)
n=100
niter <- 200  # number of iterations per simulation
rlen <- 30
# the simulations to call themselves
sims <- list(lol.sims.rtrunk, lol.sims.toep, lol.sims.rtrunk, lol.sims.fat_tails, lol.sims.qdtoep)
maxr <- c(30, 90, 30, 30, 30)
ds <- c(100, 100, 100, 1000, 100)
# additional arguments for each simulation scenario
opt_args <- list(list(), list(), list(K=3), list(rotate=TRUE), list())
sim_names = c("Trunk-2", "Toeplitz", "Trunk-3", "Fat-Tails (D=1000)", "QDA")

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

clusterExport(cl, "simulations"); clusterExport(cl, "rlen")
results <- parLapply(cl, simulations, function(sim) {
  require(lol)
  sim_dat <- do.call(sim$sim_func, sim$args)
  X <- sim_dat$X; Y <- sim_dat$Y
  results <- data.frame(sim=c(), iter=c(), alg=c(), r=c(), lhat=c())
  for (i in 1:length(algs)) {
    if (sim$sim == "QDA") {
      algs <- list(lol.project.pca, lol.project.cpca, lol.project.lrcca, lol.project.lol, lol.project.qoq)
      alg_name <- c("PCA", "cPCA", "CCA", "LOL", "QOQ")
    } else {
      algs <- list(lol.project.pca, lol.project.cpca, lol.project.lrcca, lol.project.lol)
      alg_name <- c("PCA", "cPCA", "CCA", "LOL")
    }
    rs <- round(seq(from=1, to=sim$rmax, length.out=rlen))
    for (r in rs) {
      if (alg_name %in% c("QOQ")) {
        classifier.alg=qda
      } else {
        classifier.alg=lda
      }
      tryCatch({
        xv_res <- lol.xval.eval(X, Y, alg=algs[[i]], alg.opts=list(r=r), alg.return="A", classifier=classifier.alg, k='loo')
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
resultso <- do.call(rbind, results)
results <- data.table(resultso)
# aggregate over the iterations, retaining the other factors
results.means <- aggregate(lhat ~ sim + alg + r + lhat, data = results, FUN = mean)
results_agg <- list(overall=resultso, means=results.means)
saveRDS(results_agg, 'lol_fig3_lda.rds')
stopCluster(cl)
