#' Principal Component Analysis (PCA)
#'
#' A function that performs PCA on data.
#'
#' @import irlba
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param r the rank of the projection.
#' @param ... optional args.
#' @return A list of class \code{embedding} containing the following:
#' \item{A}{\code{[d, r]} the projection matrix from \code{d} to \code{r} dimensions.}
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
#' model <- lol.project.pca(X=X, r=2)  # use pca to project into 2 dimensions
#' @export
lol.project.pca <- function(X, r, ...) {
  # mean center by the global mean
  A <- lol.utils.svd(t(Xc), r=r, center=TRUE)$A

  return(list(A=A, Xr=lol.embed(X, A)))
}

# A utility for pre-centered data to do projection faster if possible with irlba
lol.utils.svd <- function(X, r=NULL, center=FALSE, ...) {
  d <- dim(X)[2]  # dimensions of X
  if (is.null(r)) {
    r <- d
  }
  if (center) {
    X  <- sweep(X, 2, colMeans(X), '-')
  }
  X <- as.matrix(X)
  if (r < .05*d) {
    # take the svd and retain the top r left singular vectors as our components
    # using more efficient irlba if we only need a fraction of singular vecs
    svdX <- irlba::irlba(X, nv=0, nu=r)
  } else {
    svdX <- svd(X, nv=0, nu=r)
  }

  return(list(A=svdX$u, v=svdX$d))
}

#' Class PCA
#'
#' A function that performs PCA on the class-centered data. Same as low-rank LDA.
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param r the rank of the projection.
#' @param ... optional args.
#' @return A list of class \code{embedding} containing the following:
#' \item{A}{\code{[d, r]} the projection matrix from \code{d} to \code{r} dimensions.}
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
  A <- lol.utils.svd(t(Xt), r=r)$A

  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=lol.embed(X, A), cr=lol.embed(centroids, A)))
}
