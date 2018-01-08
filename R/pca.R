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
#' @return ylabs [K] vector containing the unique, ordered class labels.
#' @return centroids [K, d] centroid matrix of the unique, ordered classes.
#' @return priors [K] vector containing prior probability for the unique, ordered classes.
#' @return Xr [n, r] the data in reduced dimensionality.
#' @return cr [K, r] the centroids in reduced dimensionality.
#' @author Eric Bridgeford
#' @export
fs.project.pca <- function(X, Y, r, center=TRUE) {
  # class data
  classdat <- fs.utils.classdat(X, Y)
  priors <- classdat$priors; centroids <- classdat$centroids
  K <- classdat$K; ylabs <- classdat$ylabs
  n <- classdat$n; d <- classdat$d
  # mean center by the global mean
  if (center) {
    X <- sweep(X, 2, colMeans(X), '-')
  }

  A <- fs.utils.pca(X, r)

  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=X %*% A, cr=centroids %*% A))
}

# A utility for pre-centered data to do PCA faster.
fs.utils.pca <- function(X, r) {
  # take the svd and retain the top r left singular vectors as our components
  svd <- irlba::irlba(t(as.matrix(X)), nv=0, nu=r)
  A <- svd$u
  return(A)
}

#' Class PCA
#'
#' A function that performs PCA on the class-centered data. Same as low-rank LDA.
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
fs.project.cpca <- function(X, Y, r) {
  # class data
  classdat <- fs.utils.classdat(X, Y)
  priors <- classdat$priors; centroids <- classdat$centroids
  K <- classdat$K; ylabs <- classdat$ylabs
  n <- classdat$n; d <- classdat$d

  # subtract column means per-class
  Yidx <- sapply(Y, function(y) which(ylabs == y))
  Xt <- X - centroids[Yidx,]
  # compute the standard PCA but with the pre-centered data.
  A <- fs.utils.pca(Xt, r)

  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=X %*% A, cr=centroids %*% A))
}
