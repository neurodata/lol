## Projection Analysis
# docker pull neurodata/flashlol:0.0.3
# docker run -ti --entrypoint /bin/bash <directory/with/neurodata-public-rerf/>:/data -v <directory/lol/docs/lol-paper/data/uci>:/outputs neurodata/flashlol:0.0.3
# run as desired
#----------------------
## Load Packages
#----------------------
require(ramify)
require(parallel)
require(lolR)
require(MASS)
require(plyr)
require(dplyr)
require(randomForest)
ucipath = '/data'
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
    n <- dim(result)[1]; d <- dim(result)[2] - 1
    X <- result[, -(d+1)]; Y <- result[, d+1]
    n <- length(Y)
    k.folds <- split(1:n, 1:n.folds)
    K <- length(unique(Y))
    if (K > 10) {
      stop()
    }
    return(list(X=as.matrix(X), Y=as.vector(Y), Dataset=dset.names[i], Folds=k.folds, n=n, d=d))
  }, error = function(e) {print(i); print(e); return(NULL)})
})
# remove NULL missing elements
data <- data[which(sapply(data, function(dat) !is.null(dat)))]

experiments <- do.call(c, lapply(1:length(data), function(i) {
  n <- data[[i]]$n; d <- data[[i]]$d; k.folds <- data[[i]]$Folds
  res <- do.call(c, lapply(c("raw", "unit", "rank"), function(xfm) {
    lapply(1:n.folds, function(j) {
      return(list(Dataset.idx=i, train.idx=(1:n)[!(1:n %in% k.folds[[j]])],
                  test.idx=k.folds[[j]], rs=ceiling(seq(1, min(r.max, d), length.out=min(r.len, d)),
                  xfm=xfm, Fold=j)))
    })
  }))
  gc()
  return(res)
}))

# projection strategies of interest
proj.algs <- list(PCA=lol.project.pca, LOL=lol.project.lol, LDA=lol.project.lrlda,
                  RP=lol.project.rp, PLS=lol.project.pls, CCA=lol.project.lrcca)
# store the classifiers of interest
classifier.algs <- list(LDA=lda, RF=randomForest, RC=lol.classify.randomChance)

saveRDS(list(result=experiments), '/outputs/experiments_lol_uci.rds')
dir.create('/outputs/uci/')
#----------------------
## Algorithm Execution
#----------------------
xv.res <- mclapply(1:length(experiments), function(i) {
  exper=experiments[[i]]
  Dataset <- data[[exper$Dataset.idx]]$Dataset
  #if (file.exists(sprintf('/outputs/uci/Dataset-%s_lol_uci_fold-%s_xfm-%s.rds',
  #                        Dataset, exper$Fold, exper$xfm))) {
  #  return(readRDS(sprintf('/outputs/uci/Dataset-%s_lol_uci_fold-%s_xfm-%s.rds',
  #                         Dataset, exper$Fold, exper$xfm))$xv)
  #} else {
    tryCatch({
      X <- data[[exper$Dataset.idx]]$X; Y <- data[[exper$Dataset.idx]]$Y
      if (exper$xfm == "unit") {
        X <- scale(X, center=FALSE)
      } else if (exper$xfm == "rank") {
        X <- apply(X, 2, function(x) rank(x, ties.method="average"))
      }
      X.train <- X[exper$train.idx,,drop=FALSE]; X.test <- X[exper$test.idx,,drop=FALSE]
      Y.train <- Y[exper$train.idx]; Y.test <- Y[exper$test.idx]
      retain.cols <- apply(X.train, 2, function(x) length(unique(x))) > 1
      X.train <- X.train[,retain.cols]; X.test <- X.test[,retain.cols]
      n <- length(Y); d <- ncol(X)
      K <- length(unique(Y))
      # loop over projection strategies
      xv.res.fold <- lapply(names(proj.algs), function(proj.name) {
        proj.alg <- proj.algs[[proj.name]]
        # project using strategy proj.name
        tryCatch({
          proj.res <- do.call(proj.alg, list(X=X.train, Y=Y.train, r=min(ncol(X.train), max(exper$rs))))
          # store the training and testing projections
          op <- (list(X.train=as.matrix(proj.res$Xr), X.test=lol.embed(X.test, proj.res$A)))
        }, error=function(e) {NULL})
        lapply(exper$rs[exper$rs <= ncol(X.train)], function(r) {
          # grab the top r columns
          tryCatch({
            Xr.train <- op$X.train[,1:r, drop=FALSE]; Xr.test <- op$X.test[,1:r, drop=FALSE]
          }, error=function(e){NULL})
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
              return(data.frame(Dataset=Dataset, Algorithm=proj.name, Classifier=class.name,
                                Accuracy=acc, Misclassification=1-acc, r=r, n=n, d=d, K=K,
                                Fold=exper$Fold, xfm=exper$xfm))},
              error=function(e) {print(e); print(i); return(data.frame())})
          }) %>%
            bind_rows()
        }) %>%
          bind_rows()
      }) %>%
        bind_rows()
      rm(X.train, X.test)
      gc()
      saveRDS(list(xv=xv.res.fold, Fold=list(train=exper$train.idx, test=exper$test.idx), Y=Y),
              sprintf('/outputs/uci/Dataset-%s_lol_uci_fold-%s_xfm-%s.rds',
                      Dataset, exper$Fold, exper$xfm))
      return(xv.res.fold)
    }, error=function(e) {print(e); print(i); return(data.frame())})
  #}
}, mc.cores=detectCores()-1) %>%
  bind_rows()

saveRDS(list(result=xv.res), sprintf('/outputs/lol_uci_big.rds'))


