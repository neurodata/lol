#' Cross-Validation
#'
#' A function for performing leave-one-out cross-validation for a given model.
#' @import randomForest
#' @import MASS
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
#' @examples
#' # train model and analyze with loo validation using lda classifier
#' library(lol)
#' data <- fs.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' r=5  # embed into r=5 dimensions
#' xval.fit <- fs.xval.eval(X, Y, r, fs.project.lol, classifier='lda', k='loo')
#' @export
fs.xval.eval <- function(X, Y, r, alg, classifier='lda', k='loo') {
  Y <- factor(Y)
  n <- length(Y)
  sets <- fs.xval.split(X, Y, k=k)
  Lhat.fold <- sapply(sets, function(set) {
    mod <- do.call(alg, list(X=set$X.train, Y=set$Y.train, r=r))  # learn the projection with the algorithm specified
    X.test.proj <- set$X.test %*% mod$A  # project the data with the projection just learned
    if (classifier == 'lda') {
      liney <- MASS::lda(mod$Xr, set$Y.train)
      Yhat <- predict(liney, X.test.proj)$class
    } else if (classifier == 'rf') {
      shrubbery <- randomForest::randomForest(mod$Xr, set$Y.train)
      Yhat <- predict(shrubbery, X.test.proj)
    }
    return(1 - sum(Yhat == set$Y.test)/length(Yhat))
  })

  model <- do.call(alg, list(X=X, Y=Y, r=r))

  return(list(Lhat=mean(Lhat.fold), A=model$A, ylabs=model$ylabs, centroids=model$centroids,
              priors=model$priors, Xr=model$Xr, cr=model$cr))
}

#' Cross-Validation Data Splitter
#'
#' \code{sg.bern.xval_split_data} A function to split a dataset into
#' training and testing sets for cross validation.
#'
#' @param X [n, d] an array of input data.
#' @param Y [s] the class labels.
#' @param k='loo' the cross-validated method to perform.
#' \itemize{
#' \item{'loo'}{Leave-one-out cross validation, setting k=n.}
#' \item{isinteger(k)}{perform k-fold cross-validation with k as the number of folds.}
#' }
#' @return sets the cross-validation sets as a list, each element with an X.train, X.test, Y.train, and Y.test.
#' @author Eric Bridgeford
#' @examples
#' # prepare data for 10-fold validation
#' library(lol)
#' data <- fs.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' sets.xval.10fold <- fs.xval.split(X, Y, k=10)
#' # prepare data for loo validation
#' sets.xval.loo <- fs.xval.split(X, Y, k='loo')
#'
#' @export
fs.xval.split <- function(X, Y, k='loo') {
  Y <- factor(Y)
  n <- length(Y)
  if (k == 'loo') {
    k <- n  # loo is just xval with k=n
  }
  if (round(k) == k) {  # then xval is an integer
    samp.ids <- as.matrix(sample(1:n, n))  # the sample ids randomly permuted
    k.folds <- split(samp.ids, rep(1:k), drop=TRUE)  # split the sample ids into xval folds

    sets <- sapply(k.folds, function(fold) {
      list(X.train=X[-fold,,drop=FALSE], Y.train=Y[-fold,drop=FALSE],
           X.test=X[fold,,drop=FALSE], Y.test=Y[fold,drop=FALSE])
    }, simplify=FALSE)
  } else {
    stop("You have not entered a valid parameter for xval.")
  }
  return(sets)
}
