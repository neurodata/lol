#' Linear Optimal Low-Rank Projection (LOL)
#'
#' A function for implementing the Linear Optimal Low-Rank Projection (LOL) Algorithm.
#'
#' @param X [n, d] the data with n samples in d dimensions.
#' @param Y [n] the labels of the samples.
#' @param r the rank of the projection.
#' @return A [d, r] the projection matrix from d to r dimensions.
#' @return ylabs [C] vector containing the unique, ordered class labels.
#' @return centroids [C, d] centroid matrix of the unique, ordered classes.
#' @return priors [C] vector containing prior probability for the unique, ordered classes.
#' @author Eric Bridgeford
#' @export
fs.project.lol <- function(X, Y, r) {
  # class data
  classdat <- gs.utils.classdat(X, Y)
  priors <- classdat$priors; centroids <- classdat$centroids
  K <- classdat$C; ylabs <- classdat$ylabs
  n <- classdat$n; d <- classdat$d
  nv <- r - C

  # compute the centers per-class
  centroids <- sapply(ylabs, function(y) colMeans(X[Y==y,,drop=FALSE]))

  A <- cbind(centroids, fs.project.cpca(X))

  # orthogonalize and normalize
  A <- Matrix::qr.Q(Matrix::qr(A))
  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs))
}
