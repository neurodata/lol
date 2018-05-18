# Parallelize Stuff
#=========================#
require(MASS)
library(parallel)
require(lolR)
require(slb)
require(randomForest)
require(plyr)

no_cores = detectCores()
classifier.name <- "lda"
classifier.alg <- MASS::lda
classifier.return = 'class'
#classifier.name <- "rf"
#classifier.alg <- randomForest::randomForest
#classifier.return = NaN

rlen <- 30

# Setup Algorithms
#==========================#
algs <- list(lol.project.pca, lol.project.lrlda, lol.project.lrcca, lol.project.rp, lol.project.pls,
             lol.project.lol)
names(algs) <- c("PCA", "LRLDA", "CCA", "RP", "PLS", "LOL")
experiments <- list()
counter <- 1

data.pmlb <- slb.load.datasets(repositories="pmlb", tasks="classification", clean.invalid=TRUE, clean.ohe=10)
data.uci <- slb.load.datasets(repositories="uci", tasks="classification", clean.invalid=FALSE, clean.ohe=FALSE)
data <- c(data.uci, data.pmlb)

# Semi-Parallel
# Setup Algorithms
#=========================#

#classifier.algs <- c(lol.classify.randomGuess, MASS::lda, randomForest::randomForest)
#names(classifier.algs) <- c("RandomGuess", "LDA", "RF")
classifier.algs <- c(lol.classify.randomGuess, MASS::lda)
names(classifier.algs) <- c("RandomGuess", "LDA")
classifier.returns <- list(NULL, "class")
names(classifier.returns) <- c("RandomGuess", "LDA")

opath <- './data/'
dir.create(opath)
opath <- './data/real_data/'
dir.create(opath)
opath <- paste('./data/real_data/', classifier.name, '/', sep="")
dir.create(opath)

k = 50  # number of folds
exp <- lapply(data, function(dat) {
  tryCatch({
    if (dat$p > 100) {
      sets <- lol.xval.split(dat$X, dat$Y, k=k, rank.low=TRUE)
      return(list(sets=sets, X=dat$X, Y=dat$Y, n=dat$n, p=dat$p, K=dat$K, task=dat$task, repo=dat$repo, dataset=dat$dataset))
    } else {
      return(NULL)
    }
  }, error=function(e){return(NULL)})
})

exp <- compact(exp)

fold_rep <- data.frame(n=c(), p=c(), K=c(), task=c(), repo=c(), dataset=c(), fold=c())
for (i in 1:length(names(exp))) {
  task <- names(exp)[i]
  X <- exp[[task]]$X; Y <- exp[[task]]$Y
  n <- dim(X)[1]; d <- dim(X)[2]
  for (j in 1:(k)) {
    fold_rep <- rbind(fold_rep, data.frame(n=exp[[task]]$n, p=exp[[task]]$p, K=exp[[task]]$K, task=task,
                                    repo=exp[[task]]$repo, dataset=exp[[task]]$dataset, fold=j))
  }
}
fold_rep <- split(fold_rep, seq(nrow(fold_rep)))

results <- mclapply(fold_rep, function(fold) {
  dat <- exp[[as.character(fold$dataset)]]
  taskname <- dat$dataset
  log.seq <- function(from=0, to=15, length=rlen) {
    round(exp(seq(from=log(from), to=log(to), length.out=length)))
  }

  X <- as.matrix(dat$X); Y <- as.factor(dat$Y)
  n <- dim(X)[1]; d <- dim(X)[2]
  sets <- dat$sets
  len.set <- sapply(sets, function(set) length(set$train))
  maxr <- min(c(d - 1, min(len.set) - 1))
  sets <- list(sets[[fold$fold]])
  rs <- unique(log.seq(from=1, to=maxr, length=rlen))
  results <- data.frame(exp=c(), alg=c(), xv=c(), n=c(), ntrain=c(), d=c(), K=c(), fold=c(), r=c(), lhat=c())
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
      xv_res <- lol.xval.optimal_dimselect(X, Y, rs, algs[[i]], sets=sets,
                                           alg.opts=list(), alg.return="A", classifier=classifier.alg,
                                           classifier.return=classifier.ret, k=k)
      results <- rbind(results, data.frame(exp=taskname, alg=names(algs)[i], xv=k, n=n, ntrain=length(sets[[1]]$train), d=d, K=length(unique(Y)),
                                           fold=fold$fold, r=xv_res$folds.data$r,
                                           lhat=xv_res$folds.data$lhat, repo=dat$repo))
    }, error=function(e) {print(e); return(NULL)})
  }


  for (classifier in names(classifier.algs)) {
    for (i in 1:length(sets)) {
      tryCatch({
        set <- sets[[i]]
        model <- do.call(classifier.algs[[classifier]], list(X[set$train, ], factor(Y[set$train], levels=unique(Y[set$train]))))
        Yhat <- predict(model, X[set$test,])
        if (!is.null(classifier.returns[[classifier]])) {
          Yhat <- Yhat[[classifier.returns[[classifier]]]]
        }
        lhat <- 1 - sum(as.numeric(Yhat) == as.numeric(Y[set$test]))/length(Yhat)
        results <- rbind(results, data.frame(exp=taskname, alg=classifier, xv=k, n=n, ntrain=length(sets[[1]]$train), d=d, K=length(unique(Y)),
                                             fold=fold$fold, r=NaN, lhat=1 - max(model$priors), repo=dat$repo))
      }, error=function(e) {print(e); NULL})
    }
  }

  saveRDS(results, file=paste(opath, taskname, '_', fold$fold, '.rds', sep=""))
  return(results)
}, mc.cores=no_cores)
resultso <- do.call(rbind, results)

saveRDS(resultso, file.path(opath, paste(classifier.name, '_results.rds', sep="")))

require(stringr)

path <- './data/real_data/lda'
repo.name = 'uci'
classifier.name = 'lda'
fnames <- list.files(path, pattern='*.rds')

results <- data.frame(exp=c(), alg=c(), XV=c(), n=c(), ntrain=c(), d=c(), K=c(), fold=c(), r=c(), lhat=c(),
                       repo=c())
for (fname in fnames) {
  foldid <- strsplit(fname, '[_,.]')[[1]]
  foldid <- foldid[[length(foldid)-1]]
  dat <- readRDS(file.path(path, fname))
  dat$fold <- as.integer(foldid)
  results <- rbind(results, dat)
}
saveRDS(results, file.path(path, paste(classifier.name, '_results.rds', sep="")))
