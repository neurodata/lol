#' Data  Piling
#'
#' A function for implementing the Maximal Data Piling (MDP) Algorithm.
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param ... optional args.
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
#' \code{vignette("dp", package = "lolR")}
#'
#' @author Minh Tang and Eric Bridgeford
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.dp(X=X, Y=Y)  # use mdp to project into maximal data piling
#' @export
lol.project.dp <- function(X, Y, ...) {
  info <- lol.utils.info(X, Y)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  if (K > d) {
    stop(sprintf("The number of classes, K=%d, must be lower than the number of native dimensions, d=%d", K, d))
  }
  deltas <- lol.utils.deltas(centroids, priors)
  centroids <- t(centroids)

  # subtract column means per-class
  Yidx <- sapply(Y, function(y) which(ylabs == y))
  Xcc <- t(X - centroids[Yidx,])

  Q <- diag(d) - Xcc %*% ginv(Xcc)

  A <- qr.Q(qr(Q %*% (deltas[, 2:dim(deltas)[2], drop=FALSE])))
  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=lol.embed(X, A), cr=lol.embed(centroids, A)))
}
