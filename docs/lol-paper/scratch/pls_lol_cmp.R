# Parallelize Stuff
#=========================#
require(MASS)
library(parallel)
require(lolR)
require(slbR)
require(randomForest)
source('../plsda.R')
no_cores = detectCores() - 1
classifier.name <- "rf"
classifier.alg <- randomForest::randomForest
classifier.return = NaN

cl = makeCluster(no_cores)

# Setup Algorithms
#==========================#
algs <- list(lol.project.pls, lol.project.mpls, lol.project.opals, lol.project.qoq, lol.project.lol)
names(algs) <- c("PLS", "MPLS", "OPAL", "QOQ", "LOL")
experiments <- list()
counter <- 1

# Setup Real Data
#==========================#
rlen <- 15
ncutoff <- 1000

data <- list()

dset.names <- names(pmlb.list(task="classification")$dsets.info)
for (i in 1:length(dset.names)) {
  tryCatch({result <- pmlb.load(datasets = dset.names[i], tasks='classification', clean.nan=TRUE, clean.ohe=FALSE)
    result <- result$data[[dset.names[i]]]
    data[[dset.names[i]]] <- list(X=result$X, Y=result$Y, exp=dset.names[i])
    n <- length(result$Y)
    if (n > ncutoff) {
      k <- 10
    } else {
      k <- 'loo'
    }
    experiments[[counter]] <- list(exp=dset.names[i], k=k)
    counter <- counter + 1
  }, error = function(e) NaN)
}

# Setup Algorithms
#=========================#
opath <- './data/fig5/'
dir.create(opath)
clusterExport(cl, "data"); clusterExport(cl, "rlen")
clusterExport(cl, "experiments"); clusterExport(cl, "opath")
clusterExport(cl, "classifier.alg"); clusterExport(cl, "classifier.return")
clusterExport(cl, "classifier.name")
results <- parLapply(cl, experiments, function(exp) {
  require(lolR)
  source('../plsda.R')
  log.seq <- function(from=0, to=30, length=15) {
    round(exp(seq(from=log(from), to=log(to), length.out=length)))
  }

  algs <- list(lol.project.pls, lol.project.lol)
  alg_name <- c("PLS", "LOL")

  X <- data[[exp$exp]]$X; Y <- as.factor(data[[exp$exp]]$Y)
  n <- dim(X)[1]; d <- dim(X)[2]
  sets <- lol.xval.split(X, Y, k=exp$k)
  maxr <- min(d, 100)
  rs <- unique(log.seq(from=1, to=maxr, length=rlen))
  results <- data.frame(exp=c(), alg=c(), r=c(), n=c(), lhat=c(), fold=c())
  tryCatch({
    setTimeLimit(3000)
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
          xv_res <- lol.xval.eval(X, Y, alg=algs[[i]], sets=sets, alg.opts=list(r=r), alg.return="A", classifier=classifier.alg,
                                  classifier.return=classifier.ret, k=exp$k)
          lhat <- xv_res$Lhats
          fold <- 1:length(lhat)
          results <- rbind(results, data.frame(exp=exp$exp, alg=alg_name[i], r=r, n=n, lhat=lhat, fold=fold))
        }, error=function(e) lhat <- NaN)
      }
    }
    saveRDS(results, file=paste(opath, exp$exp, '_', classifier.name, '.rds', sep=""))
  }, error=function(e) {results <- NaN})
  return(results)
})
resultso <- do.call(rbind, results)
saveRDS(resultso, file.path(opath, paste('opal_vs_lol_', classifier.name, '.rds', sep="")))
stopCluster(cl)

# Aggregate and save
#=================================#
require(MASS)
library(parallel)
require(lolR)
require(slbR)
source('./plsda.R')

dset.names <- names(pmlb.list(task="classification")$dsets.info)
opath <- './data/fig5/'
results <- lapply(dset.names, function(dset) {
  tryCatch(
    result <- readRDS(paste(opath, dset, '.rds', sep="")), error=function(e) {return(NaN)}
  )
})
resultso <- do.call(rbind, results)
saveRDS(resultso, file.path(opath, paste('pls_v_lol_', classifier.name, '.rds', sep="")))

