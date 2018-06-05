#' Low-rank Canonical Correlation Analysis (LR-CCA)
#'
#' A function for implementing the Low-rank Canonical Correlation Analysis (LR-CCA) Algorithm.
#'
#' @import irlba
#' @param X [n, d] the data with \code{n} samples in \code{d} dimensions.
#' @param Y [n] the labels of the samples with \code{K} unique labels.
#' @param r the rank of the projection.
#' @param ... trailing args.
#' @return A list containing the following:
#' \item{\code{A}}{\code{[d, r]} the projection matrix from \code{d} to \code{r} dimensions.}
#' \item{\code{d}}{the eigen values associated with the eigendecomposition.}
#' \item{\code{ylabs}}{\code{[K]} vector containing the \code{K} unique, ordered class labels.}
#' \item{\code{centroids}}{\code{[K, d]} centroid matrix of the \code{K} unique, ordered classes in native \code{d} dimensions.}
#' \item{\code{priors}}{\code{[K]} vector containing the \code{K} prior probabilities for the unique, ordered classes.}
#' \item{\code{Xr}}{\code{[n, r]} the \code{n} data points in reduced dimensionality \code{r}.}
#' \item{\code{cr}}{\code{[K, r]} the \code{K} centroids in reduced dimensionality \code{r}.}

#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("lrcca", package = "lolR")}
#'
#' @author Eric Bridgeford and Minh Tang
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.lrcca(X=X, Y=Y, r=5)  # use lrcca to project into 5 dimensions
#' @export
lol.project.lrcca <- function(X, Y, r, ...) {
  info <- lol.utils.info(X, Y)
  priors <- info$priors; centroids <- t(info$centroids)
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }
  Yh <- lol.utils.ohe(Y)
  Xc <- X - outer(rep(1, n), colMeans(X))
  Yc <- Yh - outer(rep(1, n), colMeans(Yh))
  # covariance matrices
  S_x <- 1/n*t(Xc) %*% Xc; S_y <- 1/n*t(Yc) %*% Yc
  # inverse covariance matrices are ginverse in the low-rank case
  S_xi <- ginv(S_x); S_yi <- MASS::ginv(S_y)
  S_xy <- 1/n*t(Xc) %*% Yc
  # decompose Sxi*Sxy*Syi*Syx
  X.decomp <- lol.utils.decomp(t(S_xi %*% S_xy %*% S_yi %*% t(S_xy)), ncomp=r)

  return(list(A=X.decomp$comp, d=X.decomp$val, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=lol.embed(X, X.decomp$comp), cr=lol.embed(centroids, X.decomp$comp)))
}
