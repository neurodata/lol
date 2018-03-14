# Parallelize Stuff
#=========================#
require(lolR)
require(MASS)
library(parallel)
source('../plsda.R')
classifier.name <- "rf"
classifier.alg <- randomForest::randomForest
classifier.return = NaN

no_cores = detectCores() - 1

cl = makeCluster(no_cores)

# Setup Sims
#==========================#
require(lol)
require(MASS)

n=100
niter <- 200  # number of iterations per simulation
rlen <- 20
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
opath <- './data/fig3'
clusterExport(cl, "simulations"); clusterExport(cl, "rlen")
clusterExport(cl, "classifier.alg"); clusterExport(cl, "classifier.return")
clusterExport(cl, "classifier.name")
results <- parLapply(cl, simulations, function(sim) {
  require(lolR)
  source('./plsda.R')
  sim_dat <- do.call(sim$sim_func, sim$args)
  X <- sim_dat$X; Y <- sim_dat$Y
  results <- data.frame(sim=c(), iter=c(), alg=c(), r=c(), lhat=c())
  algs <- list(lol.project.pca, lol.project.cpca, lol.project.lrcca, lol.project.pls, lol.project.rp, lol.project.lol, lol.project.qoq)
  alg_name <- c("PCA", "LDA", "CCA", "PLS", "RP", "LOL", "QOQ")
  if (classifier.name == "lda" & sim$sim != "QDA") {
    algs <- list(lol.project.pca, lol.project.cpca, lol.project.lrcca, lol.project.pls, lol.project.rp, lol.project.lol)
    alg_name <- c("PCA", "LDA", "CCA", "RP", "PLS", "RP", "LOL")
  }

  log.seq <- function(from=0, to=30, length=15) {
    round(exp(seq(from=log(from), to=log(to), length.out=length)))
  }

  rs <- unique(round(log.seq(from=1, to=sim$rmax, length=rlen)))
  results <- data.frame(sim=c(), iter=c(), se=c(), alg=c(), r=c(), lhat=c())
  for (i in 1:length(algs)) {
    classifier.ret <- classifier.return
    if (classifier.name == "lda") {
      classifier.ret = "class"
      if (alg_name[i] == "QOQ") {
        classifier.alg=MASS::qda
        classifier.ret = "class"
      } else if (alg_name[i] == "CCA") {
        classifier.alg = lol.classify.nearestCentroid
        classifier.ret = NaN
      }
    }
    for (r in rs) {
      tryCatch({
        xv_res <- lol.xval.eval(X, Y, alg=algs[[i]], alg.opts=list(r=r), alg.return="A", classifier=classifier.alg,
                                classifier.return=classifier.ret, k='loo')
        lhat <- xv_res$Lhat
        results <- rbind(results, data.frame(sim=sim$sim, iter=sim$iter, se=var(xv_res$Lhats)/sqrt(length(Y)),
                                             alg=alg_name[i], r=r, lhat=lhat))
      }, error=function(e) lhat <- NaN)
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
results.means <- aggregate(lhat ~ sim + alg + r + lhat + se, data = results, FUN = mean)
results_agg <- list(overall=resultso, means=results.means)
saveRDS(results_agg, file.path(opath, paste('lol_fig3_', classifier.name, '.rds', sep="")))
stopCluster(cl)
