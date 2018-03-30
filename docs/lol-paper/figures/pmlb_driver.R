# Parallelize Stuff
#=========================#
require(MASS)
library(parallel)
require(lolR)
require(slbR)
require(randomForest)
no_cores = detectCores() - 1
classifier.name <- "lda"
classifier.alg <- MASS::lda
classifier.return = NaN

cl = makeCluster(no_cores)

# Setup Algorithms
#==========================#
algs <- list(lol.project.pca, lol.project.cpca, lol.project.lrcca, lol.project.rp, lol.project.pls,
             lol.project.mpls, lol.project.opal, lol.project.lol, lol.project.qoq, lol.project.plsol)
names(algs) <- c("PCA", "LDA", "CCA", "RP", "PLS", "OPAL", "MPLS", "LOL", "QOQ", "PLSOL")
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
    experiments[[counter]] <- list(exp=dset.names[i], xv=k)
    counter <- counter + 1
  }, error = function(e) NaN)
}

# Setup Algorithms
#=========================#
opath <- './data/'
dir.create(opath)
opath <- './data/real_data/'
dir.create(opath)
opath <- paste('./data/real_data/', classifier.name, '/', sep="")
dir.create(opath)
clusterExport(cl, "data"); clusterExport(cl, "rlen")
clusterExport(cl, "experiments"); clusterExport(cl, "opath")
clusterExport(cl, "classifier.alg"); clusterExport(cl, "classifier.return")
clusterExport(cl, "classifier.name"); clusterExport(cl, "algs")
results <- parLapply(cl, experiments, function(exp) {
  require(lolR)
  log.seq <- function(from=0, to=30, length=15) {
    round(exp(seq(from=log(from), to=log(to), length.out=length)))
  }

  X <- data[[exp$exp]]$X; Y <- as.factor(data[[exp$exp]]$Y)
  n <- dim(X)[1]; d <- dim(X)[2]
  maxr <- min(d, 100)
  rs <- unique(log.seq(from=1, to=maxr, length=rlen))
  sets <- lol.xval.split(X, Y, k=exp$xv)
  results <- data.frame(exp=c(), alg=c(), xv=c(), n=c(), d=c(), K=c(), fold=c(), r=c(), lhat=c())

  for (i in 1:length(algs)) {
    classifier.ret <- classifier.return
    if (classifier.name == "lda") {
      classifier.ret = "class"
      classifier.alg = MASS::lda
      if (names(algs)[i] == "QOQ") {
        classifier.alg=MASS::qda
        classifier.ret = "class"
      } else if (names(algs)[i] == "CCA") {
        classifier.alg = lol.classify.nearestCentroid
        classifier.ret = NaN
      }
    }
    tryCatch({
      xv_res <- lol.xval.optimal_r(X, Y, algs[[i]], rs, sets=sets, alg.opts=list(), alg.return="A", classifier=classifier.alg,
                                   classifier.return=classifier.ret, k=exp$xv)
      results <- rbind(results, data.frame(exp=exp$exp, alg=names(algs)[i], xv=exp$xv, n=n, d=d, K=length(unique(Y)), fold=xv_res$folds.data$fold, r=xv_res$folds.data$r,
                                           lhat=xv_res$folds.data$lhat))
    }, error=function(e) {print(e); return(NULL)})
  }
  saveRDS(results, file=paste(opath, exp$exp, '.rds', sep=""))
  return(results)
})
resultso <- do.call(rbind, results)
saveRDS(resultso, file.path(opath, paste(classifier.name, '_results.rds', sep="")))
stopCluster(cl)

# Aggregate and save
#=================================#
require(MASS)
library(parallel)
require(lolR)
require(slbR)
source('./plsda.R')

algs <- list(lol.project.pca, lol.project.cpca, lol.project.lrcca, lol.project.pls, lol.project.rp, lol.project.lol, lol.project.qoq)
names(algs) <- c("PCA", "LDA", "CCA", "PLS", "RP", "LOL", "QOQ")
experiments <- list()
counter <- 1

dset.names <- names(pmlb.list(task="classification")$dsets.info)
opath <- './data/fig5/'
results <- lapply(dset.names, function(dset) {
  tryCatch(
    result <- readRDS(paste(opath, dset, '.rds', sep="")), error=function(e) {return(NaN)}
  )
})
resultso <- do.call(rbind, results)
saveRDS(resultso, file.path(opath, paste('lol_fig5_', classifier.name, '.rds', sep="")))

