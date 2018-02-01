#' Embedding Cross Validation
#'
#' A function for performing leave-one-out cross-validation for a given embedding model.
#'
#' @importFrom MASS lda
#' @importFrom stats predict
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param alg the algorithm to use for embedding. Should be a function that accepts inputs \code{X} and \code{Y}, returning a list containing a matrix that embeds from {d} to {r < d} dimensions. Defaults to \code{lol.project.lol}.
#' @param alg.opts any extraneous options to be passed to the classifier function, as a list. Defaults to an empty list. For example, this could be the embedding dimensionality to investigate.
#' @param alg.embedding the attribute returned by \code{alg} containing the embedding matrix. Defaults to assuming that \code{alg} returns an embgedding matrix as \code{"A"}.
#' \itemize{
#' \item{!is.nan(alg.embedding)}{Assumes that \code{alg} will return a list containing an attribute, \code{alg.embedding}, a \code{[d, r]} matrix that embeds \code{[n, d]} data from \code{[d]} to \code{[r < d]} dimensions.}
#' \item{is.nan(alg.embedding)}{Assumes that \code{alg} returns a \code{[d, r]} matrix that embeds \code{[n, d]} data from \code{[d]} to \code{[r < d]} dimensions.}
#' }
#' @param classifier the classifier to use for assessing performance. The classifier should accept \code{X}, a \code{[n, d]} array as the first input, and \code{Y}, a \code{[n]} array of labels, as the first 2 arguments. The class should implement a predict function, \code{predict.classifier}, that is compatible with the \code{stats::predict} \code{S3} method. Defaults to \code{MASS::lda}.
#' @param classifier.opts any extraneous options to be passed to the classifier function, as a list. Defaults to an empty list.
#' @param classifier.return if the return type is a list, \code{class} encodes the attribute containing the prediction labels from \code{stats::predict}. Defaults to the return type of \code{MASS::lda}, \code{class}.
#' \itemize{
#' \item{!is.nan(classifier.return)}{Assumes that \code{predict.classifier} will return a list containing an attribute, \code{classifier.return}, that encodes the predicted labels.}
#' \item{is.nan(classifier.return)}{Assumes that \code{predict.classifer}} returns a \code{[n]} vector/array containing the prediction labels for \code{[n, d]} inputs.
#' }
#' @param k the cross-validated method to perform. Defaults to \code{'loo'}. See \code{\link{lol.xval.split}}
#' \itemize{
#' \item{\code{'loo'}}{Leave-one-out cross validation}
#' \item{\code{isinteger(k)}}{ perform \code{k}-fold cross-validation with \code{k} as the number of folds.}
#' }
#' @param ... trailing args.
#' @return Returns a list containing:
#' \item{Lhat}{the mean cross-validated error.}
#' \item{model}{The model returned by \code{alg} computed on all of the data.}
#' \item{classifier}{The classifier trained on all of the embedded data.}
#' \item{Lhats}{the cross-validated error for each of the \code{k}-folds.}
#' @author Eric Bridgeford
#' @examples
#' # train model and analyze with loo validation using lda classifier
#' library(lol)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' r=5  # embed into r=5 dimensions
#' # run cross-validation with the nearestCentroid method and
#' # 'eave-one-out cross-validation, which returns only
#' # prediction labels so we specify classifier.return as NaN
#' xval.fit <- lol.xval.eval(X, Y, lol.project.lol, alg.opts=list(r=r),
#'                           classifier=lol.classify.nearestCentroid,
#'                           classifier.return=NaN, k='loo')
#'
#' # train model and analyze with 5-fold validation using lda classifier
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' r=5  # embed into r=5 dimensions
#' xval.fit <- lol.xval.eval(X, Y, lol.project.lol, alg.opts=list(r=r), k=5)
#' @export
lol.xval.eval <- function(X, Y, alg, alg.opts=list(), alg.embedding="A", classifier=lda, classifier.opts=list(),
                          classifier.return="class", k='loo', ...) {
  d <- dim(X)[2]
  Y <- factor(Y)
  n <- length(Y)
  sets <- lol.xval.split(X, Y, k=k)
  Lhat.fold <- sapply(sets, function(set) {
    mod <- do.call(alg, c(list(X=set$X.train, Y=set$Y.train), alg.opts)) # learn the projection with the algorithm specified
    if (is.nan(alg.embedding)) {
      A <- mod
    } else {
      A <- mod[[alg.embedding]]
    }
    X.test.proj <- lol.embed(set$X.test, A)  # project the data with the projection just learned
    trained_classifier <- do.call(classifier, c(list(mod$Xr, set$Y.train), classifier.opts))
    if (is.nan(classifier.return)) {
      Yhat <- predict(trained_classifier, X.test.proj)
    } else {
      Yhat <- predict(trained_classifier, X.test.proj)[[classifier.return]]
    }
    return(1 - sum(Yhat == set$Y.test)/length(Yhat))
  })

  model <- do.call(alg, c(list(X=X, Y=Y), alg.opts))
  class <- do.call(classifier, c(list(model$Xr, Y), classifier.opts))

  return(list(Lhat=mean(Lhat.fold), model=model, classifier=class, Lhats=Lhat.fold))
}

#' Cross-Validation Data Splitter
#'
#' \code{sg.bern.xval_split_data} A function to split a dataset into
#' training and testing sets for cross validation.
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param k the cross-validated method to perform. Defaults to \code{'loo'}.
#' \itemize{
#' \item{\code{'loo'}}{ Leave-one-out cross validation}
#' \item{\code{isinteger(k)}}{ perform \code{k}-fold cross-validation with \code{k} as the number of folds.}
#' }
#' @param ... optional args.
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
lol.xval.split <- function(X, Y, k='loo', ...) {
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
