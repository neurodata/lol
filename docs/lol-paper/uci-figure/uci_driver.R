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
no_cores = detectCores() - 2
classifier.name <- "lda"
classifier.alg <- MASS::lda
classifier.return = 'class'
#classifier.name <- "rf"
#classifier.alg <- randomForest::randomForest
#classifier.return = NaN
ucipath = './data'

cl = makeCluster(no_cores)

# Setup Algorithms
#==========================#
algs <- list(lol.project.pca, lol.project.lrlda, lol.project.lrcca, lol.project.rp, lol.project.pls,
             lol.project.mpls, lol.project.opal, lol.project.lol, lol.project.qoq, lol.project.plsol,
             lol.project.plsolk)
names(algs) <- c("PCA", "LRLDA", "CCA", "RP", "PLS", "OPAL", "MPLS", "LOL", "QOQ", "PLSOL", "PLSOLK")
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
    k <- ceiling(n/d)
    data[[dset.names[i]]] <- list(X=X, Y=Y, exp=dset.names[i], xv=k)
  }, error = function(e) NaN)
}

classifier.algs <- list(RC=lol.classify.randomChance, LDA=MASS::lda, RF=randomForest::randomForest)

# Setup Algorithms
#=========================#
opath <- './data/'
dir.create(opath)
opath <- './data/real_data/'
dir.create(opath)
opath <- paste('./data/real_data/', classifier.name, '/', sep="")
dir.create(opath)

results <- mclapply(data, function(dat) {
  log.seq <- function(from=0, to=30, length=rlen) {
    round(exp(seq(from=log(from), to=log(to), length.out=length)))
  }

  X <- as.matrix(dat$X); Y <- as.factor(dat$Y)
  n <- dim(X)[1]; d <- dim(X)[2]
  maxr <- min(d, 100)
  rs <- unique(log.seq(from=1, to=maxr, length=rlen))
  sets <- lol.xval.split(X, Y, k=dat$xv, reverse=TRUE)
  results <- data.frame(exp=c(), alg=c(), xv=c(), n=c(), d=c(), K=c(), fold=c(), r=c(), lhat=c())

  for (i in 1:length(algs)) {
    classifier.ret <- classifier.return
    if (classifier.name == "LDA") {
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
      xv_res <- lol.xval.optimal_dimselect(X, Y, algs[[i]], rs, sets=sets, alg.opts=list(), alg.return="A", classifier=classifier.alg,
                                   classifier.return=classifier.ret, k=dat$xv)
      results <- rbind(results, data.frame(exp=dat$exp, alg=names(algs)[i], xv=dat$xv, n=n, d=d, K=length(unique(Y)),
                                           fold=xv_res$folds.data$fold, r=xv_res$folds.data$r,
                                           lhat=xv_res$folds.data$lhat))
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
        result <- rbind(result, data.frame(exp=dat$exp, alg=classifier, xv=dat$xv, n=n, d=d, K=length(unique(Y)),
                                           fold=i, r=NaN, lhat=lhat))
      }, error=function(e) {print(e); NULL})
    }
  }

  saveRDS(results, file=paste(opath, dat$exp, '.rds', sep=""))
  return(results)
}, mc.cores=no_cores)
resultso <- do.call(rbind, results)
saveRDS(resultso, file.path(opath, paste(classifier.name, '_uci_results.rds', sep="")))


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
