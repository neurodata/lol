#' Principal Component Analysis (PCA)
#'
#' A function that performs PCA on data.
#'
#' @importFrom irlba irlba
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param r the rank of the projection.
#' @param ... trailing args.
#' @return A list of class \code{embedding} containing the following:
#' \item{A}{\code{[d, r]} the projection matrix from \code{d} to \code{r} dimensions.}
#' \item{d}{\code{[r]} the signular values associated with the projection matrix \code{A}.}
#' \item{Xr}{\code{[n, r]} the \code{n} data points in reduced dimensionality \code{r}.}
#' @author Eric Bridgeford
#' @examples
#' library(lol)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.pca(X=X, r=2)  # use pca to project into 2 dimensions
#' @export
lol.project.pca <- function(X, r, ...) {
  # mean center by the global mean
  Xc  <- sweep(X, 2, colMeans(X), '-')
  svdX <- lol.utils.svd(Xc, nv=r, nu=0)

  return(list(A=svdX$v, d=svdX$d, Xr=lol.embed(X, svdX$v)))
}

#' A utility to use irlba when necessary
#' @importFrom irlba irlba
#' @author Eric Bridgeford
lol.utils.svd <- function(X, nu=0, nv=0, t=.05) {
  n <- nrow(X)
  if (nu > t*n | nv > t*n) {
    svdX <- svd(X, nu=nu, nv=nv)
  } else {
    svdX <- irlba(X, nu=nu, nv=nv)
  }
  return(svdX)
}

#' Class PCA
#'
#' A function that performs PCA on the class-centered data. Same as low-rank LDA.
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param r the rank of the projection.
#' @param ... trailing args.
#' @return A list of class \code{embedding} containing the following:
#' \item{A}{\code{[d, r]} the projection matrix from \code{d} to \code{r} dimensions.}
#' \item{d}{\code{[r]} the signular values associated with the projection matrix \code{A}.}
#' \item{ylabs}{\code{[K]} vector containing the \code{K} unique, ordered class labels.}
#' \item{centroids}{\code{[K, d]} centroid matrix of the \code{K} unique, ordered classes in native \code{d} dimensions.}
#' \item{priors}{\code{[K]} vector containing the \code{K} prior probabilities for the unique, ordered classes.}
#' \item{Xr}{\code{[n, r]} the \code{n} data points in reduced dimensionality \code{r}.}
#' \item{cr}{\code{[K, r]} the \code{K} centroids in reduced dimensionality \code{r}.}
#' @author Eric Bridgeford
#' @examples
#' library(lol)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.pca(X=X, Y=Y, r=2)  # use cpca to project into 2 dimensions
#' @export
lol.project.cpca <- function(X, Y, r, ...) {
  # class data
  classdat <- lol.utils.info(X, Y)
  priors <- classdat$priors; centroids <- t(classdat$centroids)
  K <- classdat$K; ylabs <- classdat$ylabs
  n <- classdat$n; d <- classdat$d

  # subtract column means per-class
  Yidx <- sapply(Y, function(y) which(ylabs == y))
  Xt <- X - centroids[Yidx,]
  # compute the standard projection but with the pre-centered data.
  svdX <- lol.utils.svd(Xt, nv=r, nu=0)

  return(list(A=svdX$v, d=svdX$d, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=lol.embed(X, svdX$v), cr=lol.embed(centroids, svdX$v)))
}
