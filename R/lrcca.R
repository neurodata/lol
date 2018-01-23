#' Low-rank Canonical Correlation Analysis (LR-CCA)
#'
#' A function for implementing the Low-rank Canonical Correlation Analysis (LR-CCA) Algorithm.
#'
#' @import irlba
#' @param X [n, d] the data with n samples in d dimensions.
#' @param Y [n] the labels of the samples.
#' @param r the rank of the projection.
#' @param method='full' The method to use for LR-CCA.
#' \itemize{
#' \item{'full'}{Requires O(d^2) storage, but is faster.}
#' \item{'partial'}{Requires O(n^2) storage, but is slower.}
#' }
#' @return A [d, r] the projection matrix from d to r dimensions. Only returned if method='full'.
#' @return ylabs [K] vector containing the unique, ordered class labels.
#' @return centroids [K, d] centroid matrix of the unique, ordered classes.
#' @return priors [K] vector containing prior probability for the unique, ordered classes.
#' @return Xr [n, r] the data in reduced dimensionality.
#' @return cr [K, r] the centroids in reduced dimensionality. Only returned if method='full'.
#' @author Jason Yim and Eric Bridgeford
#' @examples
#' library(lol)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.lrcca(X=X, Y=Y, r=5)  # use lrcca to project into 5 dimensions
#' @export
lol.project.lrcca <- function(X, Y, r, method='full', ...) {
  if (method == 'full') {
    outputs <- lol.project.full_lrcca(X, Y, r)
  } else if (method == 'partial') {
    outputs <- lol.project.partial_lrcca(X, Y, r)
  } else {
    stop(sprintf("You have passed an invalid method, %s. Options are 'full' and 'partial'.", method))
  }
  return(outputs)
}

lol.project.full_lrcca <- function(X, Y, r, ...) {
  info <- lol:::lol.utils.info(X, Y)
  priors <- info$priors; centroids <- t(info$centroids)
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  # hot-encoding of Y categorical variables
  Yh <- array(0, dim=c(n, K))
  # Yind is a indicator of membership in each respective class
  for (i in 1:length(ylabs)) {
    Yh[Y == ylabs[i],i] <- 1
  }

  # covariance matrices
  S_x <- cov(X); S_y <- cov(Yh)
  # inverse covariance matrices are ginverse in the low-rank case
  S_xi <- MASS::ginv(S_x); S_yi <- MASS::ginv(S_y)
  S_xy <- cov(X, y=Yh)
  # decompose Sxi*Sxy*Syi*Syx
  A <- lol.utils.pca(S_xi %*% S_xy %*% S_yi %*% t(S_xy), r)

  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=X %*% A, cr=centroids %*% A))
}

lol.project.partial_lrcca <- function(X, Y, r, ...) {
  info <- lol:::lol.utils.info(X, Y)
  priors <- info$priors; centroids <- t(info$centroids)
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d

  Y <- array(Y, dim=c(n, 1))  # force dimensionality of Y
  svd_out <- svd(X)  # svd of X
  vxvxt <- svd_out$u %*% t(svd_out$u)
  M <- svd(1/n*vxvxt %*% Y %*% t(Y))$u
  Xr <- M[, 1:r]  # dimension-reduced vectors are top r columns
  return(list(Xr=Xr, centroids=centroids, priors=priors, ylabs=ylabs))
}
