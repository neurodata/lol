#' Embedding Cross Validation
#'
#' A function for performing leave-one-out cross-validation for a given embedding model. This function produces fold-wise
#' cross-validated misclassification rates for standard embedding techniques. Users can optionally specify custom embedding techniques
#' with proper configuration of \code{alg.*} parameters and hyperparameters. Optional classifiers implementing the S3 \code{predict} function can be used
#' for classification, with hyperparameters to classifiers for determining misclassification rate specified in \code{classifier.*} parameters and hyperparameters.
#'
#' @importFrom MASS lda
#' @importFrom stats predict
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param r the number of embedding dimensions desired, where \code{r <= d}.
#' @param alg the algorithm to use for embedding. Should be a function that accepts inputs \code{X}, \code{Y}, and has a parameter for \code{alg.dimname} if \code{alg} is supervised, or just \code{X} and \code{alg.dimname} if \code{alg} is unsupervised.This algorithm should return a list containing a matrix that embeds from {d} to {r <= d} dimensions.
#' @param sets a user-defined cross-validation set. Defaults to \code{NULL}.
#' \itemize{
#' \item{\code{is.null(sets)} randomly partition the inputs \code{X} and \code{Y} into training and testing sets.}
#' \item{\code{!is.null(sets)} use a user-defined partitioning of the inputs \code{X} and \code{Y} into training and testing sets. Should be in the format of the outputs from \code{\link{lol.xval.split}}. That is, a \code{list} with each element containing \code{X.train}, an \code{[n-k][d]} subset of data to test on, \code{Y.train}, an \code{[n-k]} subset of class labels for \code{X.train}; \code{X.test}, an \code{[n-k][d]} subset of data to test the model on, \code{Y.train}, an \code{[k]} subset of class labels for \code{X.test}.}
#' }
#' @param alg.dimname the name of the parameter accepted by \code{alg} for indicating the embedding dimensionality desired. Defaults to \code{r}.
#' @param alg.opts the hyper-parameter options you want to pass into your algorithm, as a keyworded list. Defaults to \code{list()}, or no hyper-parameters.
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
#' @param k the cross-validated method to perform. Defaults to \code{'loo'}. If \code{sets} is provided, this option is ignored. See \code{\link{lol.xval.split}} for details.
#' \itemize{
#' \item{\code{'loo'} Leave-one-out cross validation}
#' \item{\code{isinteger(k)}  perform \code{k}-fold cross-validation with \code{k} as the number of folds.}
#' }
#' @param rank.low whether to force the training set to low-rank. Defaults to \code{FALSE}. If \code{sets} is provided, this option is ignored. See \code{\link{lol.xval.split}} for details.
#' \itemize{
#' \item{if \code{rank.low == FALSE}, uses default cross-validation method with standard \code{k}-fold validation. Training sets are \code{k-1} folds, and testing sets are \code{1} fold, where the fold held-out for testing is rotated to ensure no dependence of potential downstream inference in the cross-validated misclassification rates.}
#' \item{if ]code{rank.low == TRUE}, users cross-validation method with \code{ntrain = min((k-1)/k*n, d)} sample training sets, where \code{d}  is the number of dimensions in \code{X}. This ensures that the training data is always low-rank, \code{ntrain < d + 1}. Note that the resulting training sets may have \code{ntrain < (k-1)/k*n}, but the resulting testing sets will always be properly rotated \code{ntest = n/k} to ensure no dependencies in fold-wise testing.}
#' }
#' @param ... trailing args.
#' @return Returns a list containing:
#' \item{\code{lhat}}{the mean cross-validated error.}
#' \item{\code{model}}{The model returned by \code{alg} computed on all of the data.}
#' \item{\code{classifier}}{The classifier trained on all of the embedded data.}
#' \item{\code{lhats}}{the cross-validated error for each of the \code{k}-folds.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("xval", package = "lolR")}
#'
#' For  extending cross-validation techniques shown here to arbitrary embedding algorithms, see the vignette:
#' \code{vignette("extend_embedding", package = "lolR")}
#'
#' For  extending cross-validation techniques shown here to arbitrary classification algorithms, see the vignette:
#' \code{vignette("extend_classification", package = "lolR")}
#'
#' @author Eric Bridgeford
#' @examples
#' # train model and analyze with loo validation using lda classifier
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' r=5  # embed into r=5 dimensions
#' # run cross-validation with the nearestCentroid method and
#' # leave-one-out cross-validation, which returns only
#' # prediction labels so we specify classifier.return as NaN
#' xval.fit <- lol.xval.eval(X, Y, r, lol.project.lol,
#'                           classifier=lol.classify.nearestCentroid,
#'                           classifier.return=NaN, k='loo')
#'
#' # train model and analyze with 5-fold validation using lda classifier
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' xval.fit <- lol.xval.eval(X, Y, r, lol.project.lol, k=5)
#'
#' # pass in existing cross-validation sets
#' sets <- lol.xval.split(X, Y, k=2)
#' xval.fit <- lol.xval.eval(X, Y, r, lol.project.lol, sets=sets)
#' @export
lol.xval.eval <- function(X, Y, r, alg, sets=NULL, alg.dimname="r", alg.opts=list(),
                          alg.embedding="A", classifier=lda, classifier.opts=list(),
                          classifier.return="class", k='loo', rank.low=FALSE, ...) {
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
    sets <- lol.xval.split(X, Y, k=k, rank.low=rank.low)
  } else {
    lol.xval.check_xv_set(sets, n)
  }

  # hyperparameters are  the number of embedding dimensions and other options requested.
  dim.embed <- list(r)
  names(dim.embed) <- alg.dimname
  alg.hparams <- c(dim.embed, alg.opts)

  # compute fold-wise cross-validated error
  Lhat.fold <- sapply(sets, function(set) {
    mod <- do.call(alg, c(list(X=X[set$train,,drop=FALSE], Y=as.factor(Y[set$train,drop=FALSE])), alg.hparams)) # learn the projection with the algorithm specified
    if (is.nan(alg.embedding)) {
      A <- mod
    } else {
      A <- mod[[alg.embedding]]
    }
    # embed the testing points
    X.test.proj <- lol.embed(X[set$test,,drop=FALSE], A)  # project the data with the projection just learned
    # produce the desired classifier with the training data
    trained_classifier <- do.call(classifier, c(list(lol.embed(X[set$train,,drop=FALSE], A), as.factor(Y[set$train,drop=FALSE])), classifier.opts))
    # compute cross-validated error of the held-out data
    if (is.nan(classifier.return)) {
      Yhat <- predict(trained_classifier, X.test.proj)
    } else {
      Yhat <- predict(trained_classifier, X.test.proj)[[classifier.return]]
    }
    return(1 - sum(Yhat == Y[set$test,drop=FALSE])/length(Yhat))
  })

  # compute the embedding with all of the data
  model <- do.call(alg, c(list(X=X, Y=Y), alg.hparams))
  if (is.nan(alg.embedding)) {
    A <- model
  } else {
    A <- model[[alg.embedding]]
  }
  # and the classifier trained for good measure
  class <- do.call(classifier, c(list(lol.embed(X, A), Y), classifier.opts))

  return(list(lhat=mean(Lhat.fold), model=model, classifier=class, lhats=Lhat.fold))
}

nan.mean <- function(x) mean(x, na.rm=TRUE)

#' Optimal Cross-Validated Number of Embedding Dimensions
#'
#' A function for performing leave-one-out cross-validation for a given embedding model, that allows users to determine the optimal number of embedding dimensions for
#' their algorithm-of-choice. This function produces fold-wise cross-validated misclassification rates for standard embedding techniques across a specified selection of
#' embedding dimensions. Optimal embedding dimension is selected as the dimension with the lowest average misclassification rate across all folds.
#' Users can optionally specify custom embedding techniques with proper configuration of \code{alg.*} parameters and hyperparameters.
#' Optional classifiers implementing the S3 \code{predict} function can be used for classification, with hyperparameters to classifiers for
#' determining misclassification rate specified in \code{classifier.*}.
#'
#' @importFrom MASS lda
#' @importFrom stats predict
#' @importFrom stats aggregate
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels. Defaults to \code{NaN}.#' @param alg.opts any extraneous options to be passed to the classifier function, as a list. Defaults to an empty list. For example, this could be the embedding dimensionality to investigate.
#' @param rs \code{[r.n]} the embedding dimensions to investigate over, where \code{max(rs) <= d}.
#' @param alg the algorithm to use for embedding. Should be a function that accepts inputs \code{X} and \code{Y} and embedding dimension \code{r} if \code{alg} is supervised, or just \code{X} and embedding dimension \code{r} if \code{alg} is unsupervised.This algorithm should return a list containing a matrix that embeds from {d} to {r < d} dimensions.
#' @param sets a user-defined cross-validation set. Defaults to \code{NULL}.
#' \itemize{
#' \item{\code{is.null(sets)} randomly partition the inputs \code{X} and \code{Y} into training and testing sets.}
#' \item{\code{!is.null(sets)} use a user-defined partitioning of the inputs \code{X} and \code{Y} into training and testing sets. Should be in the format of the outputs from \code{\link{lol.xval.split}}. That is, a \code{list} with each element containing \code{X.train}, an \code{[n-k][d]} subset of data to test on, \code{Y.train}, an \code{[n-k]} subset of class labels for \code{X.train}; \code{X.test}, an \code{[n-k][d]} subset of data to test the model on, \code{Y.train}, an \code{[k]} subset of class labels for \code{X.test}.}
#' }
#' @param alg.dimname the name of the parameter accepted by \code{alg} for indicating the embedding dimensionality desired. Defaults to \code{r}.
#' @param alg.opts the hyper-parameter options to pass to your algorithm as a keyworded list. Defaults to \code{list()}, or no hyper-parameters. This should not include the number of embedding dimensions, \code{r}, which are passed separately in the \code{rs} vector.
#' @param alg.embedding the attribute returned by \code{alg} containing the embedding matrix. Defaults to assuming that \code{alg} returns an embgedding matrix as \code{"A"}.
#' \itemize{
#' \item \code{!is.nan(alg.embedding)} Assumes that \code{alg} will return a list containing an attribute, \code{alg.embedding}, a \code{[d, r]} matrix that embeds \code{[n, d]} data from \code{[d]} to \code{[r < d]} dimensions.
#' \item \code{is.nan(alg.embedding)} Assumes that \code{alg} returns a \code{[d, r]} matrix that embeds \code{[n, d]} data from \code{[d]} to \code{[r < d]} dimensions.
#' }
#' @param alg.structured a boolean to indicate whether the embedding matrix is structured. Provides performance increase by not having to compute the embedding matrix \code{xv} times if unnecessary. Defaults to \code{TRUE}.
#' \itemize{
#' \item \code{TRUE} assumes that if \code{Ar: R^d -> R^r} embeds from \code{d} to \code{r} dimensions and \code{Aq: R^d -> R^q} from \code{d} to \code{q > r} dimensions, that \code{Aq[, 1:r] == Ar},
#' \item \code{TRUE} assumes that if \code{Ar: R^d -> R^r} embeds from \code{d} to \code{r} dimensions and \code{Aq: R^d -> R^q} from \code{d} to \code{q > r} dimensions, that \code{Aq[, 1:r] != Ar},
#' }
#' @param classifier the classifier to use for assessing performance. The classifier should accept \code{X}, a \code{[n, d]} array as the first input, and \code{Y}, a \code{[n]} array of labels, as the first 2 arguments. The class should implement a predict function, \code{predict.classifier}, that is compatible with the \code{stats::predict} \code{S3} method. Defaults to \code{MASS::lda}.
#' @param classifier.opts any extraneous options to be passed to the classifier function, as a list. Defaults to an empty list.
#' @param classifier.return if the return type is a list, \code{class} encodes the attribute containing the prediction labels from \code{stats::predict}. Defaults to the return type of \code{MASS::lda}, \code{class}.
#' \itemize{
#' \item{\code{!is.nan(classifier.return)} Assumes that \code{predict.classifier} will return a list containing an attribute, \code{classifier.return}, that encodes the predicted labels.}
#' \item{\code{is.nan(classifier.return)} Assumes that \code{predict.classifer} returns a \code{[n]} vector/array containing the prediction labels for \code{[n, d]} inputs.}
#' }
#' @param k the cross-validated method to perform. Defaults to \code{'loo'}. If \code{sets} is provided, this option is ignored. See \code{\link{lol.xval.split}} for details.
#' \itemize{
#' \item{\code{'loo'} Leave-one-out cross validation}
#' \item{\code{isinteger(k)}  perform \code{k}-fold cross-validation with \code{k} as the number of folds.}
#' }
#' @param rank.low whether to force the training set to low-rank. Defaults to \code{FALSE}. If \code{sets} is provided, this option is ignored. See \code{\link{lol.xval.split}} for details.
#' \itemize{
#' \item{if \code{rank.low == FALSE}, uses default cross-validation method with standard \code{k}-fold validation. Training sets are \code{k-1} folds, and testing sets are \code{1} fold, where the fold held-out for testing is rotated to ensure no dependence of potential downstream inference in the cross-validated misclassification rates.}
#' \item{if ]code{rank.low == TRUE}, users cross-validation method with \code{ntrain = min((k-1)/k*n, d)} sample training sets, where \code{d}  is the number of dimensions in \code{X}. This ensures that the training data is always low-rank, \code{ntrain < d + 1}. Note that the resulting training sets may have \code{ntrain < (k-1)/k*n}, but the resulting testing sets will always be properly rotated \code{ntest = n/k} to ensure no dependencies in fold-wise testing.}
#' }
#' @param ... trailing args.
#' @return Returns a list containing:
#' \item{\code{folds.data}}{the results, as a data-frame, of the per-fold classification accuracy.}
#' \item{\code{foldmeans.data}}{the results, as a data-frame, of the average classification accuracy for each \code{r}.}
#' \item{\code{optimal.lhat}}{the classification error of the optimal \code{r}}.
#' \item{\code{optimal.r}}{the optimal number of embedding dimensions from \code{rs}}.
#' \item{\code{model}}{the model trained on all of the data at the optimal number of embedding dimensions.}
#' \item{\code{classifier}}{the classifier trained on all of the data at the optimal number of embedding dimensions.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("xval", package = "lolR")}
#'
#' For  extending cross-validation techniques shown here to arbitrary embedding algorithms, see the vignette:
#' \code{vignette("extend_embedding", package = "lolR")}
#'
#' For  extending cross-validation techniques shown here to arbitrary classification algorithms, see the vignette:
#' \code{vignette("extend_classification", package = "lolR")}
#'
#' @author Eric Bridgeford
#' @examples
#' # train model and analyze with loo validation using lda classifier
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' # run cross-validation with the nearestCentroid method and
#' # leave-one-out cross-validation, which returns only
#' # prediction labels so we specify classifier.return as NaN
#' xval.fit <- lol.xval.optimal_dimselect(X, Y, rs=c(5, 10, 15), lol.project.lol,
#'                           classifier=lol.classify.nearestCentroid,
#'                           classifier.return=NaN, k='loo')
#'
#' # train model and analyze with 5-fold validation using lda classifier
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' xval.fit <- lol.xval.optimal_dimselect(X, Y, rs=c(5, 10, 15), lol.project.lol, k=5)
#'
#' # pass in existing cross-validation sets
#' sets <- lol.xval.split(X, Y, k=2)
#' xval.fit <- lol.xval.optimal_dimselect(X, Y, rs=c(5, 10, 15), lol.project.lol, sets=sets)
#' @export
lol.xval.optimal_dimselect <- function(X, Y, rs, alg, sets=NULL, alg.dimname="r", alg.opts=list(), alg.embedding="A",
                                       alg.structured=TRUE, classifier=lda, classifier.opts=list(),
                                       classifier.return="class", k='loo', rank.low=FALSE, ...) {
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
    sets <- lol.xval.split(X, Y, k=k, rank.low=rank.low)
  } else {
    lol.xval.check_xv_set(sets, n)
  }
  # compute  the top r embedding dimensions
  max.r <- max(rs)

  # hyperparameters are  the number of embedding dimensions and other options requested.
  dim.embed <- list()
  dim.embed[[alg.dimname]] <- max.r
  alg.hparams <- c(dim.embed, alg.opts)

  Lhat.data <- lapply(1:length(sets), function(i) {
    set <- sets[[i]]
    # if the algorithm is appropriately structured, we can avoid computing on every iteration and
    # just compute on the maximum number of embedding dimensions
    if (alg.structured) {
      # learn the projection with the algorithm specified
      mod <- do.call(alg, c(list(X=X[set$train,,drop=FALSE], Y=as.factor(Y[set$train,drop=FALSE])), alg.hparams))
      # assign the embedding  dimension
      if (is.nan(alg.embedding)) {
        A <- mod
      } else {
        A <- mod[[alg.embedding]]
      }
    }
    # check fold with every desired possibility of embedding dimensions
    res.rs <- lapply(rs, function(r) {
      tryCatch({
        # if appropriately structured just take the top r dimensions of the embedding computed on max.r
        if (alg.structured) {
          A.r <- A[, 1:r]
        } else {
          # otherwise, compute A.r on the new embedding dimension every time
          alg.hparams[[alg.dimname]] <- r
          mod <- do.call(alg, c(list(X=X[set$train,,drop=FALSE], Y=as.factor(Y[set$train,drop=FALSE])), alg.hparams))
          if (is.nan(alg.embedding)) {
            A.r <- mod
          } else {
            A.r <- mod[[alg.embedding]]
          }
        }
        # embed the test points with the embedding matrix computed on the training data
        X.test.proj <- lol.embed(X[set$test,,drop=FALSE], A.r)  # project the data with the projection just learned
        # compute the trained classifier
        trained_classifier <- do.call(classifier, c(list(lol.embed(X[set$train,,drop=FALSE], A.r), as.factor(Y[set$train,drop=FALSE])), classifier.opts))
        # and compute the held-out error for particular fold
        if (is.nan(classifier.return)) {
          Yhat <- predict(trained_classifier, X.test.proj)
        } else {
          Yhat <- predict(trained_classifier, X.test.proj)[[classifier.return]]
        }
        return(data.frame(lhat=1 - sum(Yhat == Y[set$test,drop=FALSE])/length(Yhat), r=r, fold=i))
      }, error=function(e){print(e); return(NULL)})
    })
    # skip nulls
    res.rs <- res.rs[!sapply(res.rs, is.null)]
    return(do.call(rbind, res.rs))
  })

  results <- do.call(rbind, Lhat.data)

  # compute fold-wise average
  results.means <- aggregate(lhat ~ r, data = results, FUN = nan.mean)

  # find the number of embedding dimensions with lowest average error
  optimal.idx <- which(results.means$lhat == min(results.means$lhat))
  best.r <- min(results.means$r[optimal.idx])  # best is the minimum of the possible choices
  # recompute on all data with the best number of embedding dimensions
  alg.hparams[[alg.dimname]] <- best.r
  model <- do.call(alg, c(list(X=X, Y=Y), alg.hparams))
  if (is.nan(alg.embedding)) {
    A <- model
  } else {
    A <- model[[alg.embedding]]
  }
  # and train the classifier for good measure
  class <- do.call(classifier, c(list(lol.embed(X, A), Y), classifier.opts))

  return(list(folds.data=results, foldmeans.data=results.means, optimal.lhat=results.means$lhat[optimal.idx], optimal.r=best.r,
              model=model, classifier=class))
}

lol.xval.check_xv_set <- function(sets, n) {
  lapply(sets, function(set) {
    n.setmax <- max(c(set$train), c(set$test))
    if (n.setmax > n) {
      stop("You have a cross-validation set with a requested sample > total number of samples.")
    }
  })
}

#' Cross-Validation Data Splitter
#'
#' A function to split a dataset into training and testing sets for cross validation. The procedure for cross-validation
#' is to split the data into k-folds. The k-folds are then rotated individually to form a single held-out testing set the model will be validated on,
#' and the remaining (k-1) folds are used for training the developed model. Note that this cross-validation function includes functionality to be used for
#' low-rank cross-validation. In that case, instead of using the full (k-1) folds for training, we subset \code{min((k-1)/k*n, d)} samples to ensure that
#' the resulting training sets  are all low-rank. We still rotate properly over the held-out fold to ensure that the resulting testing sets
#' do not have any shared examples, which would add a complicated  dependence structure to inference we attempt to infer on the testing sets.
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param k the cross-validated method to perform. Defaults to \code{'loo'}.
#' \itemize{
#' \item{if \code{k == round(k)}, performed k-fold cross-validation.}
#' \item{if \code{k == 'loo'}, performs leave-one-out cross-validation.}
#' }
#' @param rank.low whether to force the training set to low-rank. Defaults to \code{FALSE}.
#' \itemize{
#' \item{if \code{rank == FALSE}, uses default cross-validation method with standard \code{k}-fold validation. Training sets are \code{k-1} folds, and testing sets are \code{1} fold, where the fold held-out for testing is rotated to ensure no dependence of potential downstream inference in the cross-validated misclassification rates.}
#' \item{if \code{rank == TRUE}, users cross-validation method with \code{ntrain = min((k-1)/k*n, d)} sample training sets, where \code{d}  is the number of dimensions in \code{X}. This ensures that the training data is always low-rank, \code{ntrain < d + 1}. Note that the resulting training sets may have \code{ntrain < (k-1)/k*n}, but the resulting testing sets will always be properly rotated \code{ntest = n/k} to ensure no dependencies in fold-wise testing.}
#' }
#' @param ... optional args.
#' @return sets the cross-validation sets as an object of class \code{"XV"} containing the following:
#' \item{\code{train}}{length \code{[ntrain]} vector indicating the indices of the training examples.}
#' \item{\code{test}}{ length \code{[ntest]} vector indicating the indices of the testing examples.}
#' @author Eric Bridgeford
#' @examples
#' # prepare data for 10-fold validation
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' sets.xval.10fold <- lol.xval.split(X, Y, k=10)
#'
#' # prepare data for loo validation
#' sets.xval.loo <- lol.xval.split(X, Y, k='loo')
#'
#' @export
lol.xval.split <- function(X, Y, k='loo', rank.low=FALSE, ...) {
  Y <- factor(Y)
  n <- length(Y)
  n.x <- dim(X)[1]
  d <- dim(X)[2]
  if (n != n.x) {
    stop("Your number of X samples and Y responses is not the same.")
  }
  if (k == 'loo') {
    k <- n  # loo is just xval with k=n
  }
  if (round(k) == k) {  # then xval is an integer
    samp.ids <- as.matrix(sample(1:n, n))  # the sample ids randomly permuted
    k.folds <- split(samp.ids, rep(1:k), drop=TRUE)  # split the sample ids into xval folds
    # partition X and Y appropriately into training and testing sets
    sets <- lapply(k.folds, function(fold) {
      train <- samp.ids[-fold]
      # if low-rank specified, sub-sample d elements if not already low-rank
      if (rank.low & length(train) > d) {
        train <- train[sample(1:length(train), d)]  # sample d-elements randomly
      }

      test <- samp.ids[fold]
      list(train=train, test=test)
    })
  } else {
    stop("You have not entered a valid parameter for xval.")
  }
  return(structure(sets, class="XV"))
}
