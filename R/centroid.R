#' Nearest Centroid Classifier Training
#'
#' A function that trains a classifier based on the nearest centroid.
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the \code{n} samples.
#' @param ... optional args.
#' @return A list of class \code{nearestCentroid}, with the following attributes:
#' \item{centroids}{\code{[K, d]} the centroids of each class with \code{K}  classes in \code{d} dimensions.}
#' \item{ylabs}{\code{[K]} the ylabels for each of the \code{K} unique classes, ordered.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("centroid", package = "lolR")}
#'
#' @author Eric Bridgeford
#'
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.classify.nearestCentroid(X, Y)
#' @export
lol.classify.nearestCentroid <- function(X, Y, ...) {
  # class data
  classdat <- lol.utils.info(X, Y)
  priors <- classdat$priors; centroids <- t(classdat$centroids)
  K <- classdat$K; ylabs <- classdat$ylabs
  model <-  list(centroids=centroids, ylabs=ylabs, priors=priors)
  return(structure(model, class="nearestCentroid"))
}

#' Nearest Centroid Classifier Prediction
#'
#' A function that predicts the class of points based on the nearest centroid
#' @param object An object of class \code{nearestCentroid}, with the following attributes:
#' \itemize{
#' \item{centroids}{\code{[K, d]} the centroids of each class with \code{K} classes in \code{d} dimensions.}
#' \item{ylabs}{\code{[K]} the ylabels for each of the \code{K} unique classes, ordered.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' }
#' @param X \code{[n, d]} the data to classify with \code{n} samples in \code{d} dimensions.
#' @param ... optional args.
#' @return Yhat \code{[n]} the predicted class of each of the \code{n} data point in \code{X}.
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("centroid", package = "lolR")}
#'
#' @author Eric Bridgeford
#'
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.classify.nearestCentroid(X, Y)
#' Yh <- predict(model, X)
#' @export
predict.nearestCentroid <- function(object, X, ...) {
  K <- length(object$ylabs); n <-  dim(X)[1]
  dists <- array(0, dim=c(n, K))
  for (i in 1:n) {
    for (j in 1:K) {
      dists[i, j] <- sqrt(sum((X[i,] - object$centroids[,j])^2))
    }
  }
  Yass <- apply(dists, c(1), which.min)
  Yhat <- object$ylabs[Yass]
  return(Yhat)
}
