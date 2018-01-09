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
#' @examples
#' library(fselect)
#' data <- fs.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- fs.project.pca(X=X, r=2)  # use pca to project into 2 dimensions
#' @export
fs.project.pca <- function(X, r, ...) {
  # mean center by the global mean
  Xc <- sweep(X, 2, colMeans(X), '-')
  A <- fs.utils.pca(Xc, r)

  return(list(A=A, Xr=X %*% A))
}

# A utility for pre-centered data to do PCA faster.
fs.utils.pca <- function(X, r, ...) {
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
#' @examples
#' library(fselect)
#' data <- fs.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- fs.project.pca(X=X, Y=Y, r=2)  # use cpca to project into 2 dimensions
#' @export
fs.project.cpca <- function(X, Y, r) {
  # class data
  classdat <- fselect:::fs.utils.classdat(X, Y)
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
