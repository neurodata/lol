#' Maximal Data  Piling
#'
#' A function for implementing the Maximal Data Piling (MDP) Algorithm.
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @return A list of class \code{embedding} containing the following:
#' \item{A}{\code{[d, K-1]} the projection matrix from \code{d} to \code{K-1} dimensions.}
#' \item{ylabs}{\code{[K]} vector containing the \code{K} unique, ordered class labels.}
#' \item{centroids}{\code{[K, d]} centroid matrix of the \code{K} unique, ordered classes in native \code{d} dimensions.}
#' \item{priors}{\code{[K]} vector containing the \code{K} prior probabilities for the unique, ordered classes.}
#' \item{Xr}{\code{[n, K-1]} the \code{n} data points in reduced dimensionality \code{K-1}.}
#' \item{cr}{\code{[K, K-1]} the \code{K} centroids in reduced dimensionality \code{K-1}.}
#' @author Minh Tang and Eric Bridgeford
#' @examples
#' library(lol)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.mdp(X=X, Y=Y)  # use mdp to project into maximal data piling
#' @export
lol.project.mdp <- function(X, Y, ...) {
  info <- lol:::lol.utils.info(X, Y)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  deltas <- lol:::lol.utils.deltas(centroids, priors)
  centroids <- t(centroids)

  # subtract column means per-class
  Yidx <- sapply(Y, function(y) which(ylabs == y))
  Xcc <- t(X - centroids[Yidx,])

  Q <- diag(d) - Xcc %*% ginv(Xcc)

  A <- Q %*% (deltas[, 2:dim(deltas)[2], drop=FALSE])
  return(structure(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
                        Xr=lol.embed(X, A), cr=lol.embed(centroids, A)), class="embedding"))
}
