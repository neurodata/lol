# Parallelize Stuff
#=========================#
require(MASS)
library(parallel)
require(lolR)
require(slbR)
source('./plsda.R')
no_cores = detectCores() - 1

cl = makeCluster(no_cores)

# Setup Algorithms
#==========================#
algs <- list(lol.project.pca, lol.project.cpca, lol.project.lrcca, lol.project.pls, lol.project.lol, lol.project.qoq)
names(algs) <- c("PCA", "LDA", "CCA", "PLS", "LOL", "QOQ")
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
results <- parLapply(cl, experiments, function(exp) {
  require(lolR)
  source('./plsda.R')
  log.seq <- function(from=0, to=30, length=15) {
    round(exp(seq(from=log(from), to=log(to), length.out=length)))
  }

  algs <- list(lol.project.pca, lol.project.cpca, lol.project.lrcca, lol.project.pls, lol.project.lol, lol.project.qoq)
  alg_name <- c("PCA", "LDA", "CCA", "PLS", "LOL", "QOQ")

  X <- data[[exp$exp]]$X; Y <- as.factor(data[[exp$exp]]$Y)
  n <- dim(X)[1]; d <- dim(X)[2]
  maxr <- min(d, 100)
  rs <- unique(log.seq(from=1, to=maxr, length=rlen))
  results <- data.frame(exp=c(), alg=c(), K=c(), r=c(), n=c(), lhat=c())
  tryCatch({
    setTimeLimit(1800)
    for (i in 1:length(algs)) {
      classifier.alg = MASS::lda
      classifier.return = 'class'
      if (alg_name[i] == "QOQ") {
        classifier.alg=MASS::qda
      } else if (alg_name[i] == "CCA") {
        classifier.alg = lol.classify.nearestCentroid
        classifier.return = NaN
      }
      for (r in rs) {
        tryCatch({
          xv_res <- lol.xval.eval(X, Y, alg=algs[[i]], alg.opts=list(r=r), alg.return="A", classifier=classifier.alg,
                                  classifier.return=classifier.return, k=exp$k)
          lhat <- xv_res$Lhat
          results <- rbind(results, data.frame(sim=exp$exp, se=var(xv_res$Lhats)/sqrt(length(Y)),
                                               alg=alg_name[i], r=r, lhat=lhat))
        }, error=function(e) lhat <- NaN)
      }
    }
    saveRDS(results, file=paste(opath, exp$exp, '.rds', sep=""))
  }, error=function(e) {results <- NaN})
  return(results)
})
resultso <- do.call(rbind, results)
saveRDS(resultso, 'lol_fig4pmlb_lda.rds')
stopCluster(cl)

# Aggregate and save
#=================================#
results <- lapply(dset.names, function(dset) {
  result <- readRDS(paste(opath, dset, '.rds', sep=""))
})

