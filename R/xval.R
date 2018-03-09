#' Embedding Cross Validation
#'
#' A function for performing leave-one-out cross-validation for a given embedding model.
#'
#' @importFrom MASS lda
#' @importFrom stats predict
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels. Defaults to \code{NaN}.#' @param alg.opts any extraneous options to be passed to the classifier function, as a list. Defaults to an empty list. For example, this could be the embedding dimensionality to investigate.
#' @param alg the algorithm to use for embedding. Should be a function that accepts inputs \code{X} and \code{Y} if \code{alg} is supervised, or just \code{X} if \code{alg} is unsupervised.This algorithm should return a list containing a matrix that embeds from {d} to {r < d} dimensions.
#' @param sets a user-defined cross-validation set. Defaults to \code{NULL}.
#' \itemize{
#' \item{\code{is.null(sets)} randomly partition the inputs \code{X} and \code{Y} into training and testing sets.}
#' \item{\code{!is.null(sets)} use a user-defined partitioning of the inputs \code{X} and \code{Y} into training and testing sets. Should be in the format of the outputs from \code{\link{lol.xval.split}}. That is, a \code{list} with each element containing \code{X.train}, an \code{[n-k][d]} subset of data to test on, \code{Y.train}, an \code{[n-k]} subset of class labels for \code{X.train}; \code{X.test}, an \code{[n-k][d]} subset of data to test the model on, \code{Y.train}, an \code{[k]} subset of class labels for \code{X.test}.}
#' }
#' @param alg.embedding the attribute returned by \code{alg} containing the embedding matrix. Defaults to assuming that \code{alg} returns an embgedding matrix as \code{"A"}.
#' \itemize{
#' \item \code{!is.nan(alg.embedding)} Assumes that \code{alg} will return a list containing an attribute, \code{alg.embedding}, a \code{[d, r]} matrix that embeds \code{[n, d]} data from \code{[d]} to \code{[r < d]} dimensions.
#' \item \code{is.nan(alg.embedding)} Assumes that \code{alg} returns a \code{[d, r]} matrix that embeds \code{[n, d]} data from \code{[d]} to \code{[r < d]} dimensions.
#' }
#' @param classifier the classifier to use for assessing performance. The classifier should accept \code{X}, a \code{[n, d]} array as the first input, and \code{Y}, a \code{[n]} array of labels, as the first 2 arguments. The class should implement a predict function, \code{predict.classifier}, that is compatible with the \code{stats::predict} \code{S3} method. Defaults to \code{MASS::lda}.
#' @param classifier.opts any extraneous options to be passed to the classifier function, as a list. Defaults to an empty list.
#' @param classifier.return if the return type is a list, \code{class} encodes the attribute containing the prediction labels from \code{stats::predict}. Defaults to the return type of \code{MASS::lda}, \code{class}.
#' \itemize{
#' \item{\code{!is.nan(classifier.return)} Assumes that \code{predict.classifier} will return a list containing an attribute, \code{classifier.return}, that encodes the predicted labels.}
#' \item{\code{is.nan(classifier.return)} Assumes that \code{predict.classifer} returns a \code{[n]} vector/array containing the prediction labels for \code{[n, d]} inputs.}
#' }
#' @param k the cross-validated method to perform. Defaults to \code{'loo'}. See \code{\link{lol.xval.split}}
#' \itemize{
#' \item{\code{'loo'} Leave-one-out cross validation}
#' \item{\code{isinteger(k)}  perform \code{k}-fold cross-validation with \code{k} as the number of folds.}
#' }
#' @param ... trailing args.
#' @return Returns a list containing:
#' \item{\code{Lhat}}{the mean cross-validated error.}
#' \item{\code{model}}{The model returned by \code{alg} computed on all of the data.}
#' \item{\code{classifier}}{The classifier trained on all of the embedded data.}
#' \item{\code{Lhats}}{the cross-validated error for each of the \code{k}-folds.}
#' @author Eric Bridgeford
#' @examples
#' # train model and analyze with loo validation using lda classifier
#' library(lolR)
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
#'
#' # pass in existing cross-validation sets
#' sets <- lol.xval.split(X, Y, k=2)
#' xval.fit <- lol.xval.eval(X, Y, lol.project.lol, sets=sets, alg.opts=list(r=r))
#' @export
lol.xval.eval <- function(X, Y, alg, sets=NULL, alg.opts=list(), alg.embedding="A", classifier=lda, classifier.opts=list(),
                          classifier.return="class", k='loo', ...) {
  d <- dim(X)[2]
  Y <- factor(Y)
  n <- length(Y)
  x.n <- dim(X)[1]
  if (n != x.n) {
    stop(sprintf("Your X has %d examples, but you only provide %d labels.", x.n, n))
  }
  # check that if the user specifies the cross-validation set, if so, that it is correctly set up
  # otherwise, do it for them
  if (is.null(sets)) {
    sets <- lol.xval.split(X, Y, k=k)
  } else {
    lapply(sets, function(set) {
      xtr.nk <- dim(set$X.train)[1]
      xte.k <- dim(set$X.test)[1]
      xtr.d <- dim(set$X.train)[2]
      xte.d <- dim(set$X.test)[2]
      ytr.nk <- length(set$Y.train)
      yte.k <- length(set$Y.test)
      if ((ytr.nk + yte.k != n)) {
        stop("You have a cross-validation set with ytrain.n + ytest.n != n.")
      }
      if ((xtr.nk + xte.k != n)) {
        stop("You have a cross-validation set with xtrain.n + xtest.n != n.")
      }
      if ((xtr.d != d)) {
        stop("You have a cross-validation set where xtrain.d != d.")
      }
      if ((xte.d != d)) {
        stop("You have a cross-validation set where xtest.d != d.")
      }
    })
  }

  Lhat.fold <- sapply(sets, function(set) {
    mod <- do.call(alg, c(list(X=set$X.train, Y=set$Y.train), alg.opts)) # learn the projection with the algorithm specified
    if (is.nan(alg.embedding)) {
      A <- mod
    } else {
      A <- mod[[alg.embedding]]
    }
    X.test.proj <- lol.embed(set$X.test, A)  # project the data with the projection just learned
    trained_classifier <- do.call(classifier, c(list(lol.embed(set$X.train, A), set$Y.train), classifier.opts))
    if (is.nan(classifier.return)) {
      Yhat <- predict(trained_classifier, X.test.proj)
    } else {
      Yhat <- predict(trained_classifier, X.test.proj)[[classifier.return]]
    }
    return(1 - sum(Yhat == set$Y.test)/length(Yhat))
  })

  model <- do.call(alg, c(list(X=X, Y=Y), alg.opts))
  if (is.nan(alg.embedding)) {
    A <- model
  } else {
    A <- model[[alg.embedding]]
  }
  class <- do.call(classifier, c(list(lol.embed(X, A), Y), classifier.opts))

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
#' \item \code{'loo'}  Leave-one-out cross validation
#' \item \code{isinteger(k)}  perform \code{k}-fold cross-validation with \code{k} as the number of folds.
#' }
#' @param ... optional args.
#' @return sets the cross-validation sets as a list. Each element of the list contains the following items:
#' \item{\code{X.train}}{the training data as a \code{[n - k, d]} array.}
#' \item{\code{Y.train}}{the training labels as a \code{[n - k]} vector.}
#' \item{\code{X.test}}{the testing data as a \code{[k, d]} array.}
#' \item{\code{Y.test}}{the testing labels as a \code{[k]} vector.}
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
