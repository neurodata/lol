#' Low-rank Canonical Correlation Analysis (LR-CCA)
#'
#' A function for implementing the Low-rank Canonical Correlation Analysis (LR-CCA) Algorithm.
#'
#' @import irlba
#' @param X [n, p] the data with n samples in d dimensions.
#' @param Y [n, q] the labels of the samples.
#' @param r the rank of the projection.
#' @return A [d, r] the projection matrix from d to r dimensions.
#' @return ylabs [K] vector containing the unique, ordered class labels.
#' @return centroids [K, d] centroid matrix of the unique, ordered classes.
#' @return priors [K] vector containing prior probability for the unique, ordered classes.
#' @return Xr [n, r] the data in reduced dimensionality.
#' @return cr [K, r] the centroids in reduced dimensionality.
#' @author Jason Yim and Eric Bridgeford
#' @export
fs.project.lrcca <- function(X, Y, r) {
  classdat <- fs.utils.classdat(X, Y)
  priors <- classdat$priors; centroids <- classdat$centroids
  K <- classdat$K; ylabs <- classdat$ylabs
  n <- classdat$n; d <- classdat$d

  Yind <- array(0, dim=c(n, K))
  # Yind is a indicator of membership in each respective class
  for (i in 1:length(ylabs)) {
    Yind[Y == y,i] <- 1
  }
  X <- as.array(X)
  # covariance matrices
  S_x <- cov(X); S_y <- cov(Yind)
  # inverse covariance matrices are ginverse in the low-rank case
  S_xi <- MASS::ginv(S_x); S_yi <- MASS::ginv(S_y)
  S_xy <- cov(X, y=Yind); S_yx <- cov(Yind, y=X)
  # decompose
  A <- fs.utils.pca(S_xi %*% S_xy %*% S_yi %*% S_yx, r)

  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=X %*% A, cr=centroids %*% A))
}
