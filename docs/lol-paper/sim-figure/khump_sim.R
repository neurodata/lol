# Parallelize Stuff
#=========================#
require(lolR)
require(MASS)
library(parallel)
#classifier.name <- "rf"
#classifier.alg <- randomForest::randomForest
#classifier.return = NaN
classifier.name <- 'lda'
classifier.alg <- MASS::lda
classifier.return <- "class"

no_cores = detectCores() - 1

n=128
niter <- 50  # number of iterations per simulation
rlen <- 20
# the simulations to call themselves
sims <- list(lol.sims.khump)#, lol.sims.kident)
maxr <- c(40)

# additional arguments for each simulation scenario
opt_args <- list(list())# ,list())
sim_names = c("Hump-K")#, "Identity-K")

algs <- list(lol.project.pca, lol.project.lrlda, lol.project.lrcca, lol.project.rp, lol.project.pls,
             lol.project.lol)
names(algs) <- c("PCA", "LRLDA", "CCA", "RP", "PLS", "LOL")
alg.opts=list(list(), list(), list(), list(), list(), list())
names(alg.opts) <- c("PCA", "LRLDA", "CCA", "RP", "PLS", "LOL")

simulations <- list()
counter <- 1
d=100
for (i in 1:length(sims)) {
  for (j in 1:niter) {
    for (k in seq(2, 10, length.out=5)) {
      simulations[[counter]] <- list(sim_func=sims[[i]], args=c(list(n, d, K=k)),
                                     rmax=maxr, sim=sim_names[i], iter=j)
      counter <- counter + 1
    }
  }
}

#=========================#
opath <- './data/sims'
results <- mclapply(simulations, function(sim) {

  sim_dat <- do.call(sim$sim_func, sim$args)
  X <- sim_dat$X; Y <- sim_dat$Y
  results <- data.frame(sim=c(), iter=c(), alg=c(), r=c(), lhat=c())

  log.seq <- function(from=0, to=30, length=rlen) {
    round(exp(seq(from=log(from), to=log(to), length.out=length)))
  }

  rs <- unique(round(log.seq(from=1, to=sim$rmax, length=rlen)))
  result <- do.call(rbind, lapply(1:length(algs), function(i) {
    classifier.ret <- classifier.return
    if (classifier.name == "lda") {
      classifier.ret = "class"
      classifier.alg = MASS::lda
    }
    if (sim$sim == "Cross") {  # for cross simulation, use QDA
      classifier.alg=MASS::qda
      classifier.ret = "class"
    }
    if (names(algs)[i] == "CCA") {  # CCA produces turrible embeddings that dont work w LDA as they are singular
      classifier.alg = lol.classify.nearestCentroid
      classifier.ret = NaN
    }
    tryCatch({
      xv_res <- lol.xval.optimal_dimselect(X, Y, rs, algs[[names(algs)[i]]], alg.opts=alg.opts[[names(algs)[i]]],
                                           alg.return="A", classifier=classifier.alg,
                                           classifier.return=classifier.ret, k='loo')
      data.frame(sim=sim$sim, iter=sim$iter, fold=xv_res$folds.data$fold, alg=names(algs)[i], K=sim$args$K,
                                           r=xv_res$folds.data$r, lhat=xv_res$folds.data$lhat)
    }, error=function(e) data.frame(sim=sim$sim, iter=sim$iter, fold=xv_res$folds.data$fold, alg=names(algs)[i], K=sim$args$K,
                                    r=xv_res$folds.data$r, lhat=NaN))
  }))
  return(result)
}, mc.cores=no_cores)

# Aggregate and save
#=================================#
resultso <- do.call(rbind, results)
saveRDS(resultso, file.path(opath, paste('lol_sims_khump_', classifier.name, '.rds', sep="")))

