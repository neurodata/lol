#' Embedding Cross Validation
#'
#' A function for performing leave-one-out cross-validation for a given embedding model.
#'
#' @import randomForest
#' @import MASS
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param r the number of dimensions to project the data onto. Should have \code{r < d}.
#' @param alg=lol.project.lol the algorithm to use for embedding. Should be a function returning something of class \code{embedding}.
#' @param classifier='lda' the classifier to use for assessing performance.
#' \itemize{
#' \item{'lda'}{ Use the lda classifier for assessing performance. \code{\link[MASS]{lda}}}
#' \item{'rf'}{ Use the random forest classifier for assessing performance. \code{\link[randomForest]{randomForest}}}
#' \item{'cent'}{ Use the nearest centroid classifier for assessing performance \code{\link{lol.classifier.nearestCentroid}}}
#' }
#' @param k='loo' the cross-validated method to perform. \code{\link{lol.xval.split}}
#' \itemize{
#' \item{\code{'loo'}}{Leave-one-out cross validation}
#' \item{\code{isinteger(k)}}{ perform \code{k}-fold cross-validation with \code{k} as the number of folds.}
#' }
#' @return Returns a list containing:
#' \item{Lhat}{the mean cross-validated error.}
#' \item{A}{\code{[d, r]} the projection matrix from \code{d} to \code{r} dimensions.}
#' \item{ylabs}{\code{[K]} vector containing the \code{K} unique, ordered class labels.}
#' \item{centroids}{\code{[K, d]} centroid matrix of \code{K} the unique, ordered classes in \code{d} dimensions.}
#' \item{priors}{\code{[K]} vector containing the \code{K} prior probabilities for the unique, ordered classes.}
#' \item{Xr}{\code{[n, r]} the \code{n} data points  in reduced dimensionality \code{r}.}
#' \item{cr}{\code{[K, r]} the \code{K} centroids in reduced dimensionality \code{r}.}
#' \item{Lhats}{the cross-validated error for each of the \code{k}-folds.}
#' @author Eric Bridgeford
#' @examples
#' # train model and analyze with loo validation using lda classifier
#' library(lol)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' r=5  # embed into r=5 dimensions
#' xval.fit <- lol.xval.eval(X, Y, r, lol.project.lol, classifier='lda', k='loo')
#'
#' # train model and analyze with loo val5-fold validation using rf classifier
#' library(lol)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' r=5  # embed into r=5 dimensions
#' xval.fit <- lol.xval.eval(X, Y, r, lol.project.lol, classifier='rf', k='5')
#' @export
lol.xval.eval <- function(X, Y, r, alg, classifier='lda', k='loo', ...) {
  d <- dim(X)[2]
  if (r > d) {
    stop(sprintf("You have specified to reduce dimensionality to %d, but native dimensionality is %d.", r, d))
  }
  Y <- factor(Y)
  n <- length(Y)
  sets <- lol.xval.split(X, Y, k=k)
  Lhat.fold <- sapply(sets, function(set) {
    mod <- do.call(alg, list(X=set$X.train, Y=set$Y.train, r=r)) # learn the projection with the algorithm specified
    X.test.proj <- set$X.test %*% mod$A  # project the data with the projection just learned
    if (classifier == 'lda') {
      liney <- MASS::lda(mod$Xr, set$Y.train)
      Yhat <- predict(liney, X.test.proj)$class
    } else if (classifier == 'rf') {
      shrubbery <- randomForest::randomForest(mod$Xr, set$Y.train)
      Yhat <- predict(shrubbery, X.test.proj)
    } else if (classifier == 'cent') {
      droid <- lol::lol.classify.nearestCentroid(mod$Xr, set$Y.train)
      Yhat <- predict(droid, X.test.proj)
    }
    return(1 - sum(Yhat == set$Y.test)/length(Yhat))
  })

  model <- do.call(alg, list(X=X, Y=Y, r=r))

  return(list(Lhat=mean(Lhat.fold), A=model$A, ylabs=model$ylabs, centroids=model$centroids,
              priors=model$priors, Xr=model$Xr, cr=model$cr, Lhats=Lhat.fold))
}

#' Cross-Validation Data Splitter
#'
#' \code{sg.bern.xval_split_data} A function to split a dataset into
#' training and testing sets for cross validation.
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param k='loo' the cross-validated method to perform.
#' \itemize{
#' \item{\code{'loo'}}{ Leave-one-out cross validation}
#' \item{\code{isinteger(k)}}{ perform \code{k}-fold cross-validation with \code{k} as the number of folds.}
#' }
#' @return sets the cross-validation sets as a list, each element with an \code{X.train}, \code{X.test}, \code{Y.train}, and \code{Y.test}.
#' @author Eric Bridgeford
#' @examples
#' # prepare data for 10-fold validation
#' library(lol)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' sets.xval.10fold <- lol.xval.split(X, Y, k=10)
#'
#' # prepare data for loo validation
#' sets.xval.loo <- lol.xval.split(X, Y, k='loo')
#'
#' @export
lol.xval.split <- function(X, Y, k='loo') {
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
