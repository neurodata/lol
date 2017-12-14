#' Cross-Validation
#'
#' A function for performing leave-one-out cross-validation for a given model.
#' @param X [n, d] the data as n samples in d dimensions.
#' @param Y [n] the labels for each for each of the n samples.
#' @param r the number of dimensions to project the data onto.
#' @param alg=fs.project.lol the algorithm to use for feature selection.
#' @param classifier='lda' the classifier to use for assessing performance.
#' \itemize{
#' \item{'lda'}{Use the lda classifier for assessing performance.}
#' \item{'rf'}{Use the random forest classifier for assessing performance.}
#' }
#' @param k='loo' the cross-validated method to perform.
#' \itemize{
#' \item{'loo'}{Leave-one-out cross validation}
#' \item{isinteger(k)}{perform k-fold cross-validation with k as the number of folds.}
#' }
#' @return Lhat the cross-validated error.
#' @return A [d, r] the projection matrix.
#' @return ylabs [K] vector containing the unique, ordered class labels.
#' @return centroids [K, d] centroid matrix of the unique, ordered classes.
#' @return priors [K] vector containing prior probability for the unique, ordered classes.
#' @return Xr [n, r] the data in reduced dimensionality.
#' @return cr [K, r] the centroids in reduced dimensionality.
#' @author Eric Bridgeford
#' @export
fs.eval.xval <- function(X, Y, r, alg, classifier='lda', k='loo') {
  Y <- factor(Y)
  n <- length(Y)
  if (k == 'loo') {
    k <- n  # loo is just xval with k=n
  }
  if (round(k) == k) {  # then xval is an integer
    samp.ids <- as.matrix(sample(1:n, n))  # the sample ids randomly permuted
    k.folds <- split(samp.ids, rep(1:k), drop=TRUE)  # split the sample ids into xval folds
    Lhat.fold <- sapply(k.folds, function(fold) {
      X.train <- X[-fold,,drop=FALSE]
      Y.train <- Y[-fold,drop=FALSE]
      X.test <- X[fold,,drop=FALSE]
      Y.test <- Y[fold,drop=FALSE]
      mod <- do.call(alg, list(X.train, Y.train, r))  # learn the projection with the algorithm specified
      X.test.proj <- X.test %*% mod$A  # project the data with the projection just learned
      if (classifier == 'lda') {
        Yhat <- fs.lda(mod$Xr, Y.train, X.test.proj)
      } else if (classifier == 'rf') {
        shrubbery <-randomForest(mod$Xr, Y.train)
        Yhat <- predict(shrubbery, X.test.proj)
      }
      return(sum(Yhat == Y.test)/length(fold))
    })
  } else {
    stop("You have not entered a valid parameter for xval.")
  }
  model <- do.call(alg, list(X, Y, r))

  return(list(Lhat=1 - mean(Lhat.fold), A=model$A, ylabs=model$ylabs, centroids=model$centroids,
              priors=model$priors, Xr=model$Xr, cr=model$cr))
}
