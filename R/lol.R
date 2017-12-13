#' Linear Optimal Low-Rank Projection (LOL)
#'
#' A function for implementing the Linear Optimal Low-Rank Projection (LOL) Algorithm.
#'
#' @param X [n, d] the data with n samples in d dimensions.
#' @param Y [n] the labels of the samples.
#' @param r the rank of the projection.
#' @return A [d, r] the projection matrix from d to r dimensions.
#' @return ylabs [K] vector containing the unique, ordered class labels.
#' @return centroids [K, d] centroid matrix of the unique, ordered classes.
#' @return priors [K] vector containing prior probability for the unique, ordered classes.
#' @return Xr [n, r] the data in reduced dimensionality.
#' @return cr [K, r] the centroids in reduced dimensionality.
#' @author Eric Bridgeford
#' @export
fs.project.lol <- function(X, Y, r) {
  # class data
  classdat <- fs.utils.classdat(X, Y)
  priors <- classdat$priors; centroids <- classdat$centroids
  K <- classdat$K; ylabs <- classdat$ylabs
  n <- classdat$n; d <- classdat$d
  nv <- r - K
  if (nv > 0) {
    A <- cbind(t(centroids), fs.project.cpca(X, Y, nv)$A)
  } else {
    A <- t(centroids)
  }

  # orthogonalize and normalize
  A <- Matrix::qr.Q(Matrix::qr(A))
  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=X %*% A, cr=centroids %*% A))
}
