#' Randomly Guessing Classifier Training
#'
#' A function that predicts by randomly guessing based on the pmf of the class priors. Functionality consistent
#' with the standard R prediction interface so that one can compute the "guess" accuracy
#' with minimal modification of other classification scripts.
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the \code{n} samples.
#' @param ... optional args.
#' @return A list of class \code{randomGuess}, with the following attributes:
#' \item{ylabs}{\code{[K]} the ylabels for each of the \code{K} unique classes, ordered.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' @author Eric Bridgeford
#'
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.classify.randomGuess(X, Y)
#' @export
lol.classify.randomGuess <- function(X, Y, ...) {
  model <- lol.classify.rand(X, Y)
  return(structure(model, class="randomGuess"))
}

#' Randomly Guessing Classifier Prediction
#'
#' A function that predicts by randomly guessing based on the pmf of the class priors. Functionality consistent
#' with the standard R prediction interface so that one can compute the "guess" accuracy
#' with minimal modification of other classification scripts.
#' @param object An object of class \code{randomGuess}, with the following attributes:
#' \itemize{
#' \item{ylabs}{\code{[K]} the ylabels for each of the \code{K} unique classes, ordered.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' }
#' @param X \code{[n, d]} the data to classify with \code{n} samples in \code{d} dimensions.
#' @param ... optional args.
#' @return Yhat \code{[n]} the predicted class of each of the \code{n} data point in \code{X}.
#' @author Eric Bridgeford
#'
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.classify.randomGuess(X, Y)
#' Yh <- predict(model, X)
#' @export
predict.randomGuess <- function(object, X, ...) {
  K <- length(object$ylabs); n <-  dim(X)[1]
  Yass <- sample(x=seq(1:K), prob=object$priors, size=n, replace=TRUE)
  Yhat <- object$ylabs[Yass]
  return(Yhat)
}

#' Randomly Chance Classifier Training
#'
#' A function that predicts the maximally present class in the dataset. Functionality consistent
#' with the standard R prediction interface so that one can compute the "chance" accuracy
#' with minimal modification of other classification scripts.
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the \code{n} samples.
#' @param ... optional args.
#' @return A list of class \code{randomGuess}, with the following attributes:
#' \item{ylabs}{\code{[K]} the ylabels for each of the \code{K} unique classes, ordered.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' @author Eric Bridgeford
#'
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.classify.randomChance(X, Y)
#' @export
lol.classify.randomChance <- function(X, Y, ...) {
  model <- lol.classify.rand(X, Y)
  return(structure(model, class="randomChance"))
}

#' Randomly Chance Classifier Prediction
#'
#' A function that predicts the maximally present class in the dataset. Functionality consistent
#' with the standard R prediction interface so that one can compute the "chance" accuracy
#' with minimal modification of other classification scripts.
#' @param object An object of class \code{randomChance}, with the following attributes:
#' \itemize{
#' \item{ylabs}{\code{[K]} the ylabels for each of the \code{K} unique classes, ordered.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' }
#' @param X \code{[n, d]} the data to classify with \code{n} samples in \code{d} dimensions.
#' @param ... optional args.
#' @return Yhat \code{[n]} the predicted class of each of the \code{n} data point in \code{X}.
#' @author Eric Bridgeford
#'
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.classify.randomChance(X, Y)
#' Yh <- predict(model, X)
#' @export
predict.randomChance <- function(object, X, ...) {
  K <- length(object$ylabs); n <-  dim(X)[1]
  return(rep(object$ylabs[which.max(object$priors)], n))
}

#' Random Classifier Utility
#'
#' A function for random classifiers.
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the \code{n} samples.
#' @param ... optional args.
#' @return A structure, with the following attributes:
#' \item{ylabs}{\code{[K]} the ylabels for each of the \code{K} unique classes, ordered.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' @author Eric Bridgeford
lol.classify.rand <- function(X, Y, ...) {
  # class data
  classdat <- lol.utils.info(X, Y)
  priors <- classdat$priors; centroids <- t(classdat$centroids)
  K <- classdat$K; ylabs <- classdat$ylabs
  model <-  list(ylabs=ylabs, priors=priors)
  return(model)
}
