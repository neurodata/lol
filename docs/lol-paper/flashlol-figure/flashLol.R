#' Principal Component Analysis (PCA)
#'
#' A function that performs PCA on data.
#'
#' @importFrom irlba irlba
#' @importFrom robustbase colMedians
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param r the rank of the projection.
#' @param xfm whether to transform the variables before taking the SVD.
#' \itemize{
#' \item{FALSE}{apply no transform to the variables.}
#' \item{'unit'}{unit transform the variables, defaulting to centering and scaling to mean 0, variance 1. See \link[base]{scale} for details and optional arguments to be passed with \code{xfm.opts}.}
#' \item{'log'}{log-transform the variables, for use-cases such as having high variance in larger values. Defaults to natural logarithm. See \link[base]{log} for details and optional arguments to be passed with \code{xfm.opts}.}
#' \item{'rank'}{rank-transform the variables. Defalts to breaking ties with the average rank of the tied values. See \link[base]{rank} for details and optional arguments to be passed with \code{xfm.opts}.}
#' \item{c(opt1, opt2, etc.)}{apply the transform specified in opt1, followed by opt2, etc.}
#' }
#' @param xfm.opts optional arguments to pass to the \code{xfm} option specified. Should be a numbered list of lists, where \code{xfm.opts[[i]]} corresponds to the optional arguments for \code{xfm[i]}. Defaults to the default options for each transform scheme.
#' @param robust whether to perform PCA on a robust estimate of the covariance matrix or not. Defaults to \code{FALSE}.
#' @param ... trailing args.
#' @return A list containing the following:
#' \item{\code{Xr}}{\code{[n, r]} the \code{n} data points in reduced dimensionality \code{r}.}
#' \item{\code{d}}{the eigen values associated with the eigendecomposition.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("pca", package = "lolR")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.pca(X=X, r=2)  # use pca to project into 2 dimensions
#' @export
flashx.pca <- function(X, r, ...) {
  # mean center by the column mean
  d <- dim(X)[2]
  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }
  # center the data
  Xc  <- sweep(X, 2, colMeans(X), '-')
  X.decomp <- flashx.decomp(Xc, ncomp=r)

  return(list(Xr = flashx.embed(X, X.decomp$comp), d=X.decomp$val, A=X.decomp$comp))
}

flashx.embed <- function(X, A) {
  return(as.matrix(X %*% A))
}

#' A utility to use irlba when necessary
#' @importFrom irlba irlba
#' @importFrom robust covRob
#' @param X the data to compute the svd of.
#' @param ncomp the number of left singular vectors to retain.
#' @param t the threshold of percent of singular vals/vecs to use irlba.
#' @param xfm whether to transform the variables before taking the SVD.
#' \itemize{
#' \item{FALSE}{apply no transform to the variables.}
#' \item{'unit'}{unit transform the variables, defaulting to centering and scaling to mean 0, variance 1. See \link[base]{scale} for details and optional args.}
#' \item{'log'}{log-transform the variables, for use-cases such as having high variance in larger values. Defaults to natural logarithm. See \link[base]{log} for details and optional args.}
#' \item{'rank'}{rank-transform the variables. Defalts to breaking ties with the average rank of the tied values. See \link[base]{rank} for details and optional args.}
#' \item{c(opt1, opt2, etc.)}{apply the transform specified in opt1, followed by opt2, etc.}
#' }
#' @param xfm.opts optional arguments to pass to the \code{xfm} option specified. Should be a numbered list of lists, where \code{xfm.opts[[i]]} corresponds to the optional arguments for \code{xfm[i]}. Defaults to the default options for each transform scheme.
#' @param robust whether to use a robust estimate of the covariance matrix when taking PCA. Defaults to \code{FALSE}.
#' @return the svd of X.
#' @author Eric Bridgeford
flashx.decomp <- function(X, ncomp=0) {
  svdX <- fm.svd(X, nu=0, nv=ncomp)
  decomp=list(comp=svdX$v, val=svdX$d)
  return(decomp)
}

#' Low-Rank Linear Discriminant Analysis (LRLDA)
#'
#' A function that performs LRLDA on the class-centered data. Same as class-conditional PCA.
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param r the rank of the projection.
#' @param xfm whether to transform the variables before taking the SVD.
#' \itemize{
#' \item{FALSE}{apply no transform to the variables.}
#' \item{'unit'}{unit transform the variables, defaulting to centering and scaling to mean 0, variance 1. See \link[base]{scale} for details and optional args.}
#' \item{'log'}{log-transform the variables, for use-cases such as having high variance in larger values. Defaults to natural logarithm. See \link[base]{log} for details and optional args.}
#' \item{'rank'}{rank-transform the variables. Defalts to breaking ties with the average rank of the tied values. See \link[base]{rank} for details and optional args.}
#' \item{c(opt1, opt2, etc.)}{apply the transform specified in opt1, followed by opt2, etc.}
#' }
#' @param xfm.opts optional arguments to pass to the \code{xfm} option specified. Should be a numbered list of lists, where \code{xfm.opts[[i]]} corresponds to the optional arguments for \code{xfm[i]}. Defaults to the default options for each transform scheme.
#' @param robust whether to use a robust estimate of the covariance matrix when taking PCA. Defaults to \code{FALSE}.
#' @param ... trailing args.
#' @return A list containing the following:
#' \item{\code{d}}{the eigen values associated with the eigendecomposition.}
#' \item{\code{Xr}}{\code{[n, r]} the \code{n} data points in reduced dimensionality \code{r}.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("lrlda", package = "lolR")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.lrlda(X=X, Y=Y, r=2)  # use cpca to project into 2 dimensions
#' @export
flashx.lrlda <- function(X, Y, r, ...) {
  # class data
  n <- length(Y); d <- ncol(X)
  classdat <- Y %>%
    table()
  K = length(classdat); ylabs <- unique(Y)

  centroids <- fm.mapply.col(
    # compute the group-wise sums
    fm.groupby(X, 2, fm.as.factor(Y, K), fm.bo.add),
    # normalize each column by the class mean for that column
    classdat, fm.bo.div)

  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }

  # subtract column means per-class
  Yidx <- sapply(Y, function(y) which(ylabs == y))
  # form class-conditional data matrix
  Xc <- X - fm.get.rows(centroids, Yidx)
  # compute the standard projection but with the pre-centered data
  X.decomp <- flashx.decomp(Xc, ncomp=r)

  return(list(Xr = flashx.embed(X, X.decomp$comp), d=X.decomp$val, A=X.decomp$comp))
}

flashx.deltas <- function(centroids, priors) {
  d <- nrow(centroids); K <- length(priors)
  # compute the rank-K difference space as deltas(i) = mus(i) - mus(0) where the mus are ordered
  # by decreasing prior
  deltas <- array(0, dim=c(d, K))
  srt_prior <- sort(priors, decreasing=TRUE, index.return=TRUE)$ix
  gr_mix <- srt_prior[1]
  deltas[[1]] <- fm.get.cols(centroids, str_prior[1])
  for (i in 2:K) {
    deltas[[i]] <- fm.get.cols(centroids, str_prior[i]) - deltas[[1]]
  }
  deltas <- fm.cbind.list(deltas)
  return(deltas)
}

flashx.mdiff <- function(X, Y, ...) {
  # class data
  n <- length(Y); d <- ncol(X)
  classdat <- Y %>%
    table()
  K = length(classdat); ylabs <- unique(Y)
  priors <- classdat/n

  centroids <- fm.mapply.col(
    # compute the group-wise sums
    fm.groupby(X, 2, fm.as.factor(Y, K), fm.bo.add),
    # normalize each column by the class mean for that column
    classdat, fm.bo.div)

  deltas <- flashx.deltas(centroids, priors)
  return(list(Xr=flashx.embed(X, deltas), A=A))
}

flashx.rp <- function(X, r, ...) {
  n <- nrow(X); d <- ncol(X)
  A <- 1/r*fm.rnorm.matrix(d, r, 0, 1, in.mem=TRUE)
  return(list(Xr=flashx.embed(X, A), A=A))
}

flashx.lrcca <- function(X, Y, r, ...) {
  n <- lrngth(Y); d=ncol(X)
  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }
  Yh <- lolR::lol.utils.ohe(Y)
  S_y <- cov(Yh)
  S_yinv <- MASS::ginv(S_y)
  # covariance matrices cov(X)
  #Xc <- X - outer(rep(1, n), colMeans(X))
  if (nrow(X) < ncol(X)) {
    Xc <- sweep(X, 2, colMeans(X), "-")
    #S_x <- 1/(n-1)*t(Xc) %*% Xc;
    X_cov_mul <- function(vec, extra) 1/(n-1) * t(Xc) %*% (Xc %*% vec)
    res <- fm.eigen(X_cov_mul, k=nrow(X), n=ncol(X), which="LM", sym=TRUE)
  } else {
    S_x <- cov(X)
    res <- eigen(S_x)
  }
  sigma <- res$values[res$values > 1e-6]
  U <- res$vectors[,1:length(sigma)]

  # inverse covariance matrices are ginverse in the low-rank case
  # S_xi <- U %*% diag(1/sigma) %*% t(U)
  # S_xi <- MASS::ginv(S_x);
  # S_xy <- 1/(n-1)*t(Xc) %*% Yc
  S_xy <- cov(X, fm.as.matrix(Yh))

  # decompose Sxi*Sxy*Syi*Syx
  # A <- lol.utils.pca(S_xi %*% S_xy %*% S_yi %*% t(S_xy), r, trans=FALSE)
  Z <- (t(U) %*% S_xy) %*% (S_yi %*% t(S_xy))
  Z <- t(Z)
  U <- U %*% diag(1/sigma)
  A <- tcrossprod_pca(U, Z, r)

  return(list(Xr=lol.embed(X, A), A=A))
}

tcrossprod_pca <- function(V, U, r) {
  V <- fm.as.matrix(V)
  U <- fm.as.matrix(U)
  mul <- function(vec, extra) {
    vec <- as.vector(U %*% (t(V) %*% vec))
    vec <- V %*% (t(U) %*% vec)
  }
  res <- fm.eigen(mul, k=r, n=nrow(V), which="LM", sym=TRUE)
  res$vectors
}
