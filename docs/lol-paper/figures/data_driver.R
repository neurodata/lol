# Parallelize Stuff
#=========================#
require(MASS)
library(parallel)
require(lolR)
require(slb)
require(randomForest)

no_cores = detectCores() - 5
classifier.name <- "lda"
classifier.alg <- MASS::lda
classifier.return = 'class'
#classifier.name <- "rf"
#classifier.alg <- randomForest::randomForest
#classifier.return = NaN
ucipath = './data'

rlen <- 30

# Setup Algorithms
#==========================#
algs <- list(lol.project.pca, lol.project.lrlda, lol.project.lrcca, lol.project.rp, lol.project.pls,
             lol.project.mpls, lol.project.opal, lol.project.lol, lol.project.qoq, lol.project.plsol,
             lol.project.plsolk)
names(algs) <- c("PCA", "LRLDA", "CCA", "RP", "PLS", "OPAL", "MPLS", "LOL", "QOQ", "PLSOL", "PLSOLK")
experiments <- list()
counter <- 1

data.pmlb <- slb.load.datasets(repositories="pmlb", tasks="classification", clean.invalid=TRUE, clean.ohe=10)
data.uci <- slb.load.datasets(repositories="uci", tasks="classification", clean.invalid=FALSE, clean.ohe=FALSE)
data <- c(data.pmlb, data.uci)

# Semi-Parallel
# Setup Algorithms
#=========================#

classifier.algs <- c(lol.classify.randomGuess, MASS::lda, randomForest::randomForest)
names(classifier.algs) <- c("RandomGuess", "LDA", "RF")

opath <- './data/'
dir.create(opath)
opath <- './data/real_data/'
dir.create(opath)
opath <- paste('./data/real_data/', classifier.name, '/', sep="")
dir.create(opath)

cl = makeCluster(no_cores)
clusterExport(cl, "data"); clusterExport(cl, "rlen")
clusterExport(cl, "opath")
clusterExport(cl, "classifier.alg"); clusterExport(cl, "classifier.return")
clusterExport(cl, "classifier.name"); clusterExport(cl, "algs")
clusterExport(cl, "classifier.algs")
results <- parLapply(cl, data, function(dat) {
  require(lolR)
  taskname <- dat$dataset
  log.seq <- function(from=0, to=30, length=rlen) {
    round(exp(seq(from=log(from), to=log(to), length.out=length)))
  }

  X <- as.matrix(dat$X); Y <- as.factor(dat$Y)
  n <- dim(X)[1]; d <- dim(X)[2]
  if (d > 50) {
    # if the problem is not full-rank, make it full-rank by doing clever cross-validation
    if (d < n) {
      k = ceiling(n/d)
      if (k == n/d) {
        k <- k + 1
      }
      sets <- lol.xval.split(X, Y, k=k, reverse=TRUE)
      if (k < 20) {
        k = 20
      }
      sets <- sets[names(sets)[1:k]]
    } else {
      k=20
      sets <- lol.xval.split(X, Y, k=k)
    }
    len.set <- sapply(sets, function(set) length(set$train))
    maxr <- min(c(d - 1, min(len.set) - 1, 100))
    rs <- unique(log.seq(from=1, to=maxr, length=rlen))
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
        xv_res <- lol.xval.optimal_dimselect(X, Y, rs, algs[[i]], sets=sets, alg.opts=list(), alg.return="A", classifier=classifier.alg,
                                             classifier.return=classifier.ret, k=k)
        results <- rbind(results, data.frame(exp=taskname, alg=names(algs)[i], xv=k, n=n, d=d, K=length(unique(Y)),
                                             fold=xv_res$folds.data$fold, r=xv_res$folds.data$r,
                                             lhat=xv_res$folds.data$lhat, repo=dat$repo))
      }, error=function(e) {print(e); return(NULL)})
    }


    for (classifier in names(classifier.algs)) {
      for (i in 1:length(sets)) {
        tryCatch({
          set <- sets[[i]]
          model <- do.call(classifier.algs[[classifier]], list(set$X.train, factor(set$Y.train, levels=unique(set$Y.train))))
          Yhat <- predict(model, set$X.test)
          if (!is.null(classifier.return[[classifier]])) {
            Yhat <- Yhat[[classifier.return[[classifier]]]]
          }
          lhat <- 1 - sum(as.numeric(Yhat) == as.numeric(set$Y.test))/length(Yhat)
          result <- rbind(result, data.frame(exp=taskname, alg=classifier, xv=k, n=n, d=d, K=length(unique(Y)),
                                             fold=i, r=NaN, lhat=lhat, repo=dat$repo))
        }, error=function(e) {print(e); NULL})
      }
    }

    saveRDS(results, file=paste(opath, taskname, '.rds', sep=""))
    return(results)
  } else {
    return(NULL)
  }
})
resultso <- do.call(rbind, results)
saveRDS(resultso, file.path(opath, paste(classifier.name, '_results.rds', sep="")))
stopCluster(cl)


## Fully-Parallel
# Parallelize Stuff
#=========================#
require(MASS)
library(parallel)
require(lolR)
require(slb)
require(randomForest)

no_cores = detectCores() - 5
classifier.name <- "lda"
classifier.alg <- MASS::lda
classifier.return = 'class'
#classifier.name <- "rf"
#classifier.alg <- randomForest::randomForest
#classifier.return = NaN
ucipath = './data'

rlen <- 30

# Setup Algorithms
#==========================#
algs <- list(lol.project.pca, lol.project.lrlda, lol.project.lrcca, lol.project.rp, lol.project.pls,
             lol.project.mpls, lol.project.opal, lol.project.lol, lol.project.qoq, lol.project.plsol,
             lol.project.plsolk)
names(algs) <- c("PCA", "LRLDA", "CCA", "RP", "PLS", "OPAL", "MPLS", "LOL", "QOQ", "PLSOL", "PLSOLK")
experiments <- list()
counter <- 1

data.pmlb <- slb.load.datasets(repositories="pmlb", tasks="classification", clean.invalid=TRUE, clean.ohe=10)
data.uci <- slb.load.datasets(repositories="uci", tasks="classification", clean.invalid=FALSE, clean.ohe=FALSE)
data <- c(data.pmlb, data.uci)

# Semi-Parallel
# Setup Algorithms
#=========================#

classifier.algs <- c(lol.classify.randomGuess, MASS::lda, randomForest::randomForest)
names(classifier.algs) <- c("RandomGuess", "LDA", "RF")

opath <- './data/'
dir.create(opath)
opath <- './data/real_data/'
dir.create(opath)
opath <- paste('./data/real_data/', classifier.name, '/', sep="")
dir.create(opath)
