require(slb)
n.folds <- 50

saveRDS(lapply(unique(as.character(readRDS('../data/real_data/lda_results.rds')$exp)), function(dset) {
  tryCatch({
    result <- slb.load.datasets(datasets = dset, tasks='classification', clean.nan=TRUE, clean.ohe=FALSE)
    Y <- result[[dset]]$Y
    folds <- split(1:length(Y), rep(1:n.folds), drop=TRUE)  # split the sample ids into xval folds
    lapply(1:length(folds), function(i) {
      fold <- folds[[i]]
      Y.train <- Y[-fold]; Y.test <- Y[fold]
      y.guess <- names(sort(table(Y.train),decreasing=TRUE)[1])
      Accuracy <- mean(Y.test == y.guess)
      return(data.frame(Dataset=dset, Classifier="RC", Accuracy=Accuracy,
                      Er.Rt=1-Accuracy, Fold=i, n=length(Y)))
    }) %>%
      bind_rows()
  }, error=function(e) {return(data.frame())})
}) %>%
  bind_rows(), '../data/real_data/chance.rds')
