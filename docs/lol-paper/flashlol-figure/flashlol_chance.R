saveRDS(lapply(list.files('../data/flashlol/Datasets/'), function(fname) {
  result <- readRDS(sprintf('../data/flashlol/Datasets/%s', fname))
  Dataset <- strsplit(fname, '-|_')[[1]][2]
  folds <- result$k.folds; Y <- result$Y
  lapply(1:length(folds), function(i) {
    fold <- folds[[i]]
    Y.train <- Y[-fold]; Y.test <- Y[fold]
    y.guess <- names(sort(table(Y.train),decreasing=TRUE)[1])
    Accuracy <- mean(Y.test == y.guess)
    return(data.frame(Dataset=Dataset, Classifier="RC", Accuracy=Accuracy,
                      Er.Rt=1-Accuracy, Fold=i, n=length(Y)))
  }) %>%
    bind_rows()
}) %>%
  bind_rows(), '../data/flashlol/chance.rds')
