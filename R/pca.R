#' Principal Component Analysis (PCA)
#'
#' A function that performs PCA on data.
#'
#' @import irlba
#' @param X [n, d] the data with n samples in d dimensions.
#' @param Y [n] the labels of the samples.
#' @param r the rank of the projection.
#' @param center=TRUE whether to center the data before applying PCA.
#' @return A [d, r] the projection matrix from d to r dimensions.
#' @return ylabs [C] vector containing the unique, ordered class labels.
#' @return centroids [C, d] centroid matrix of the unique, ordered classes.
#' @return priors [C] vector containing prior probability for the unique, ordered classes.
#' @author Eric Bridgeford
#' @export
fs.project.pca <- function(Xt, Y, r, center=TRUE) {
  # class data
  classdat <- gs.utils.classdat(X, Y)
  priors <- classdat$priors; centroids <- classdat$centroids
  K <- classdat$C; ylabs <- classdat$ylabs
  n <- classdat$n; d <- classdat$d
  # mean center by the global mean
  if (center) {
    X <- X - colMeans(X)
  }
  # take the svd and retain the top r left singular vectors as our components
  svd <- irlba::irlba(t(as.matrix(X)), nv=0, nu=r)
  A <- svd$u

  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs))
}

#' Class PCA
#'
#' A function that performs PCA on the class-centered data.
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
fs.project.cpca <- function(X, Y, r) {
  # class data
  classdat <- gs.utils.classdat(X, Y)
  priors <- classdat$priors; centroids <- classdat$centroids
  K <- classdat$C; ylabs <- classdat$ylabs
  n <- classdat$n; d <- classdat$d

  # subtract column means per-class
  Xt <- sapply(ylabs, function(y) X[Y==y,,drop=FALSE] - centroids[ylabs == y,,drop=FALSE])
  # compute the standard PCA but with the pre-centered data.
  A <- fs.project.pca(Xt, Y, r, center=FALSE)

  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs))
}
