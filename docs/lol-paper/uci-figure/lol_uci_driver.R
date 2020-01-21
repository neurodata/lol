## Projection Analysis

#----------------------
## Load Packages
#----------------------
require(ramify)
require(parallel)
require(FlashR)
require(MASS)
require(dplyr)
require(plyr)
require(randomForest)
ucipath = '../data'
n.folds <- 50
r.len <- 20
r.max <- 100

#----------------------
## Prepare Data
#----------------------
dset.names <- readLines(file.path(ucipath, 'processed/names.txt'))
data <- lapply(1:length(dset.names), function(i) {
  tryCatch({
    result <- read.csv(file.path(ucipath, 'processed/data', paste(dset.names[i], '.csv', sep="")), header=FALSE)
    n <- dim(result)[1]; d <- dim(result)[2]
    X <- result[, -d]; Y <- result[, d]
    n <- length(Y)
    k.folds <- split(1:n, 1:n.folds)
    return(list(X=X, Y=Y, Dataset=dset.names[i], Folds=k.folds, n=n, d=d))
  }, error = function(e) return(NULL))
})
# remove NULL missing elements
data <- data[which(sapply(data, is.null))]

runner.dat <- do.call(rbind, lapply(1:length(data), function(i) {
  n <- data[[i]]$n; d <- data[[i]]$d; k.folds <- data[[i]]$Folds
  lapply(1:n.folds, function(j) {
    return(list(Dataset.idx=i, train.idx=(1:n)[!(1:n %in% k.folds[[j]])],
                test.idx=k.folds[[j]], rs=seq(1, min(r.max, d), length.out=min(r.len, d))))
  })
  gc()
}))

# projection strategies of interest
proj.algs <- list(PCA=lol.project.pca, LOL=lol.project.lol, LDA=lol.project.lrlda,
                  RP=lol.project.rp, PLS=lol.project.pls, CCA=lol.project.cca)
# store the classifiers of interest
classifier.algs <- list(LDA=lda, RF=randomForest, RC=lol.classify.randomChance)
# compute which ids will be part of which folds
k.folds <- split(1:length(Y), rep(1:n.folds), drop=TRUE)  # split the sample ids into xval folds

#----------------------
## Algorithm Execution
#----------------------
xv.res <- mclapply(runner.dat, function(exp) {
  X <- data[[exp$Dataset.idx]]$X; Y <- data[[exp$Dataset.idx]]$Y
  X.train <- X[exp$train.idx,,drop=FALSE]; X.test <- X[exp$test.idx,,drop=FALSE]
  Y.train <- Y[exp$train.idx]; Y.test <- Y[exp$test.idx]
  n <- length(Y); d <- ncol(X)
  # loop over projection strategies
  xv.res.fold <- lapply(names(proj.algs), function(proj.name) {
    proj.alg <- proj.algs[[proj.name]]
    # project using strategy proj.name
    proj.res <- do.call(proj.alg, list(X=X.train, Y=Y.train, r=max(exp$rs)))
    # store the training and testing projections
    op <- (list(X.train=as.matrix(proj.res$Xr), X.test=flashx.embed(X.test, proj.res$A)))
    lapply(exp$rs, function(r) {
      # grab the top r columns
      Xr.train <- op$X.train[,1:r, drop=FALSE]; Xr.test <- op$X.test[,1:r, drop=FALSE]
      lapply(names(classifier.algs), function(class.name) {
        tryCatch({
          class.alg <- classifier.algs[[class.name]]
          # train classifier of interest
          trained.classifier <- do.call(class.alg, list(Xr.train, as.factor(Y.train)))
          Y.hat <- predict(trained.classifier, Xr.test)  # make predictions
          if (class.name == "LDA") {
            Y.hat <- Y.hat$class  # MASS::lda puts the class labels in this attribute
          }
          # compute accuracy
          acc <- mean(Y.hat == Y.test)
          # return as a data frame
          return(data.frame(Dataset=data[[exp$Dataset.idx]]$Dataset, Algorithm=proj.name, Classifier=class.name,
                            Accuracy=acc, Misclassification=1-acc, r=r, n=n, d=d, Fold=j))},
          error=function(e) {print(e); return(NULL)})
      }) %>%
        bind_rows()
    }) %>%
      bind_rows()
  }) %>%
    bind_rows()
  rm(X.train, X.test)
  gc()
  saveRDS(list(xv=xv.res.fold), sprintf('/data/uci/Dataset-%s_lol_uci_fold-%s.rds', dataset, j))
  return(xv.res.fold)
}, mc.cores=detectCores()-1) %>%
  bind_rows()

saveRDS(list(result=xv.res, k.folds=k.folds, Y=Y), sprintf('/data/Dataset-%s_lol_uci.rds', dataset))
