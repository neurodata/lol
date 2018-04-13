function terminal() {
  sudo mkdir /data/
  sudo chmod -R 777 /data/
  aws s3 cp --no-sign-request s3://neurodata-public-rerf/uci/processed.zip ./data
  cd ./data/
  unzip processed.zip
}

# Parallelize Stuff
#=========================#
require(MASS)
library(parallel)
require(lolR)
require(slbR)
require(randomForest)
no_cores = detectCores() - 5
classifier.name <- "lda"
classifier.alg <- MASS::lda
classifier.name <- "rf"
classifier.alg <- randomForest::randomForest
classifier.return = NaN
ucipath = './data'

cl = makeCluster(no_cores)

# Setup Algorithms
#==========================#
algs <- list(lol.project.pca, lol.project.cpca, lol.project.lrcca, lol.project.rp, lol.project.pls,
             lol.project.mpls, lol.project.opal, lol.project.lol, lol.project.qoq, lol.project.plsol,
             lol.project.plsolk)
names(algs) <- c("PCA", "LDA", "CCA", "RP", "PLS", "OPAL", "MPLS", "LOL", "QOQ", "PLSOL", "PLSOLK")
experiments <- list()
counter <- 1

# Setup Real Data
#==========================#
rlen <- 20
ncutoff <- 1

data <- list()

dset.names <- readLines(file.path(ucipath, 'processed/names.txt'))
for (i in 1:length(dset.names)) {
  tryCatch({
    result <- read.csv(file.path(ucipath, 'processed/data', paste(dset.names[i], '.csv', sep="")), header=FALSE)
    n <- dim(result)[1]; d <- dim(result)[2]
    X <- result[, -d]; Y <- result[, d]
    n <- length(Y)
    if (n > ncutoff) {
      k <- 20
    } else {
      k <- 'loo'
    }
    data[[dset.names[i]]] <- list(X=X, Y=Y, exp=dset.names[i], xv=k)
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
results <- parLapply(cl, data, function(dat) {
  require(lolR)
  log.seq <- function(from=0, to=30, length=rlen) {
    round(exp(seq(from=log(from), to=log(to), length.out=length)))
  }

  X <- as.matrix(dat$X); Y <- as.factor(dat$Y)
  n <- dim(X)[1]; d <- dim(X)[2]
  maxr <- min(d, 100)
  rs <- unique(log.seq(from=1, to=maxr, length=rlen))
  sets <- lol.xval.split(X, Y, k=dat$xv)
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
                                   classifier.return=classifier.ret, k=dat$xv)
      results <- rbind(results, data.frame(exp=dat$exp, alg=names(algs)[i], xv=dat$xv, n=n, d=d, K=length(unique(Y)), fold=xv_res$folds.data$fold, r=xv_res$folds.data$r,
                                           lhat=xv_res$folds.data$lhat))
    }, error=function(e) {print(e); return(NULL)})
  }

  for (i in 1:length(sets)) {
    set <- sets[[i]]
    model <- lol.classify.randomGuess(set$X.train, set$Y.train)
    Yhat <- predict(model, set$X.test)
    lhat <- 1 - sum(Yhat == set$Y.test)/length(Yhat)
    results <- rbind(results, data.frame(exp=dat$exp, alg='RandomGuess', xv=dat$xv, n=n, d=d, K=length(unique(Y)), fold=i, r=NaN, lhat=lhat))
  }

  saveRDS(results, file=paste(opath, dat$exp, '.rds', sep=""))
  return(results)
})
resultso <- do.call(rbind, results)
saveRDS(resultso, file.path(opath, paste(classifier.name, '_uci_results.rds', sep="")))
stopCluster(cl)


classifier.name <- "lda"
opath <- paste('./data/real_data/', classifier.name, '/', sep="")
exp_names = list.files(opath)
results <- lapply(exp_names, function(dset) {
  tryCatch(
    result <- readRDS(paste(opath, dset, sep="")), error=function(e) {return(NaN)}
  )
})
resultso <- do.call(rbind, results)
saveRDS(resultso, file.path(opath, paste(classifier.name, '_uci_results.rds', sep="")))
