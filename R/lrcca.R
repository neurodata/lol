#' Low-rank Canonical Correlation Analysis (LR-CCA)
#'
#' A function for implementing the Low-rank Canonical Correlation Analysis (LR-CCA) Algorithm.
#'
#' @import irlba
#' @param X [n, d] the data with \code{n} samples in \code{d} dimensions.
#' @param Y [n] the labels of the samples with \code{K} unique labels.
#' @param r the rank of the projection.
#' @param ... trailing args.
#' @return A list of containing the following:
#' \item{A}{\code{[d, r]} the projection matrix from \code{d} to \code{r} dimensions.}
#' \item{d}{\code{[r]} the signular values associated with the projection matrix \code{A}.}
#' \item{ylabs}{\code{[K]} vector containing the \code{K} unique, ordered class labels.}
#' \item{centroids}{\code{[K, d]} centroid matrix of the \code{K} unique, ordered classes in native \code{d} dimensions.}
#' \item{priors}{\code{[K]} vector containing the \code{K} prior probabilities for the unique, ordered classes.}
#' \item{Xr}{\code{[n, r]} the \code{n} data points in reduced dimensionality \code{r}.}
#' \item{cr}{\code{[K, r]} the \code{K} centroids in reduced dimensionality \code{r}.}
#' @author Eric Bridgeford and Minh Tang
#' @examples
#' library(lol)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.lrcca(X=X, Y=Y, r=5)  # use lrcca to project into 5 dimensions
#' @export
lol.project.lrcca <- function(X, Y, r, ...) {
  info <- lol.utils.info(X, Y)
  priors <- info$priors; centroids <- t(info$centroids)
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  # hot-encoding of Y categorical variables
  Yh <- array(0, dim=c(n, K))
  # Yind is a indicator of membership in each respective class
  for (i in 1:length(ylabs)) {
    Yh[Y == ylabs[i],i] <- 1
  }
  Xc <- X - outer(rep(1, n), colMeans(X))
  Yc <- Yh - outer(rep(1, n), colMeans(Yh))
  # covariance matrices
  S_x <- 1/n*t(Xc) %*% Xc; S_y <- 1/n*t(Yc) %*% Yc
  # inverse covariance matrices are ginverse in the low-rank case
  S_xi <- ginv(S_x); S_yi <- MASS::ginv(S_y)
  S_xy <- 1/n*t(Xc) %*% Yc
  # decompose Sxi*Sxy*Syi*Syx
  svdX <- lol.utils.svd(S_xi %*% S_xy %*% S_yi %*% t(S_xy), nu=r)

  return(list(A=svdX$u, d=svdX$d, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=lol.embed(X, svdX$u), cr=lol.embed(centroids, svdX$u)))
}

lol.project.full_lrcca <- function(X, Y, r, ...) {
  info <- lol.utils.info(X, Y)
  priors <- info$priors; centroids <- t(info$centroids)
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  # hot-encoding of Y categorical variables
  Yh <- array(0, dim=c(n, K))
  # Yind is a indicator of membership in each respective class
  for (i in 1:length(ylabs)) {
    Yh[Y == ylabs[i],i] <- 1
  }
  Xc <- X - outer(rep(1, n), colMeans(X))
  Yc <- Yh - outer(rep(1, n), colMeans(Yh))
  # covariance matrices
  S_x <- 1/n*t(Xc) %*% Xc; S_y <- 1/n*t(Yc) %*% Yc
  # inverse covariance matrices are ginverse in the low-rank case
  S_xi <- ginv(S_x); S_yi <- MASS::ginv(S_y)
  S_xy <- 1/n*t(Xc) %*% Yc
  # decompose Sxi*Sxy*Syi*Syx
  svdX <- lol.utils.svd(S_xi %*% S_xy %*% S_yi %*% t(S_xy), nu=r)

  return(list(A=svdX$u, d=svdX$d, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=lol.embed(X, svdX$u), cr=lol.embed(centroids, svdX$u)))
}

lol.project.partial_lrcca <- function(X, Y, r, ...) {
  info <- lol.utils.info(X, Y)
  priors <- info$priors; centroids <- t(info$centroids)
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d

  # hot-encoding of Y categorical variables
  Yh <- array(0, dim=c(n, K))
  # Yind is a indicator of membership in each respective class
  for (i in 1:length(ylabs)) {
    Yh[Y == ylabs[i],i] <- 1
  }
  Xc <- X - outer(rep(1, n), colMeans(X))
  Yc <- Yh - outer(rep(1, n), colMeans(Yh))

  svd_out <- svd(Xc, nv=0)  # svd of X
  vxvxt <- svd_out$u %*% t(svd_out$u)
  Xr <- svd(1/n*vxvxt %*% Yc %*% t(Yc), nv=0, nu=r)$u
  return(list(Xr=Xr, centroids=centroids, priors=priors, ylabs=ylabs))
}
