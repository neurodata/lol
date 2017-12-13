#' Low-rank Canonical Correlation Analysis (LR-CCA)
#'
#' A function for implementing the Low-rank Canonical Correlation Analysis (LR-CCA) Algorithm.
#'
#' @param X [n, p] the data with n samples in d dimensions.
#' @param Y [n, q] the labels of the samples.
#' @param r the rank of the projection.
#' @return A [d, r] the projection matrix from d to r dimensions.
#' @return ylabs [K] vector containing the unique, ordered class labels.
#' @return centroids [K, d] centroid matrix of the unique, ordered classes.
#' @return priors [K] vector containing prior probability for the unique, ordered classes.
#' @return Xr [n, r] the data in reduced dimensionality.
#' @return cr [K, r] the centroids in reduced dimensionality.
#' @author Jason Yim
#' @export
fs.project.lrcca <- function(X, Y, r) {
  # class data
  classdat <- fs.utils.classdat(X, Y)
  priors <- classdat$priors; centroids <- classdat$centroids
  K <- classdat$K; ylabs <- classdat$ylabs
  n <- classdat$n; d <- classdat$d

  # canonical correlation
  cxy <- stats::cancor(X, Y)
  A <-cxy$xcoef[,1:r]
  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=X %*% A, cr=centroids %*% A))
}
