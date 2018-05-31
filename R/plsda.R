#' Partial Least-Squares (PLS)
#'
#' A function for implementing the Partial Least-Squares (PLS) Algorithm.
#'
#' @importFrom pls plsr
#' @param X [n, d] the data with \code{n} samples in \code{d} dimensions.
#' @param Y [n] the labels of the samples with \code{K} unique labels.
#' @param r the rank of the projection.
#' @param ... trailing args.
#' @return A list containing the following:
#' \item{\code{A}}{\code{[d, r]} the projection matrix from \code{d} to \code{r} dimensions.}
#' \item{\code{ylabs}}{\code{[K]} vector containing the \code{K} unique, ordered class labels.}
#' \item{\code{centroids}}{\code{[K, d]} centroid matrix of the \code{K} unique, ordered classes in native \code{d} dimensions.}
#' \item{\code{priors}}{\code{[K]} vector containing the \code{K} prior probabilities for the unique, ordered classes.}
#' \item{\code{Xr}}{\code{[n, r]} the \code{n} data points in reduced dimensionality \code{r}.}
#' \item{\code{cr}}{\code{[K, r]} the \code{K} centroids in reduced dimensionality \code{r}.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("pls", package = "lolR")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.pls(X=X, Y=Y, r=5)  # use pls to project into 5 dimensions
#' @export
lol.project.pls <- function(X, Y, r, ...) {
  info <- lol.utils.info(X, Y)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d

  Yh <- lol.utils.ohe(Y)
  tmpData <- data.frame(n=paste("row", 1:nrow(Yh), sep=""))
  tmpData$Y <- Yh
  tmpData$X <- X

  A <- plsr(Y ~ X, data=tmpData, ncomp=r)$projection
  centroids <- t(centroids)
  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=lol.embed(X, A), cr=lol.embed(centroids, A)))
}
