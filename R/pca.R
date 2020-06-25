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
#' \item{'unit'}{unit transform the variables, defaulting to centering and scaling to mean 0, variance 1. See \code{\link[base]{scale}} for details and optional arguments to be passed with \code{xfm.opts}.}
#' \item{'log'}{log-transform the variables, for use-cases such as having high variance in larger values. Defaults to natural logarithm. See \code{\link[base:Log]{log}} for details and optional arguments to be passed with \code{xfm.opts}.}
#' \item{'rank'}{rank-transform the variables. Defalts to breaking ties with the average rank of the tied values. See \code{\link[base]{rank}} for details and optional arguments to be passed with \code{xfm.opts}.}
#' \item{c(opt1, opt2, etc.)}{apply the transform specified in opt1, followed by opt2, etc.}
#' }
#' @param xfm.opts optional arguments to pass to the \code{xfm} option specified. Should be a numbered list of lists, where \code{xfm.opts[[i]]} corresponds to the optional arguments for \code{xfm[i]}. Defaults to the default options for each transform scheme.
#' @param robust whether to perform PCA on a robust estimate of the covariance matrix or not. Defaults to \code{FALSE}.
#' @param ... trailing args.
#' @return A list containing the following:
#' \item{\code{A}}{\code{[d, r]} the projection matrix from \code{d} to \code{r} dimensions.}
#' \item{\code{d}}{the eigen values associated with the eigendecomposition.}
#' \item{\code{Xr}}{\code{[n, r]} the \code{n} data points in reduced dimensionality \code{r}.}
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
lol.project.pca <- function(X, r, xfm=FALSE, xfm.opts=list(), robust=FALSE,...) {
  # mean center by the column mean
  d <- dim(X)[2]
  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }
  # subtract means
  if (!robust) {
    Xc  <- sweep(X, 2, colMeans(X), '-')
  } else {
    Xc  <- sweep(X, 2, colMedians(X), '-')
  }
  X.decomp <- lol.utils.decomp(Xc, xfm=xfm, xfm.opts=xfm.opts, ncomp=r, robust=robust)

  return(list(A=X.decomp$comp, d=X.decomp$val, Xr=lol.embed(X, X.decomp$comp)))
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
#' \item{'unit'}{unit transform the variables, defaulting to centering and scaling to mean 0, variance 1. See \code{\link[base]{scale}} for details and optional args.}
#' \item{'log'}{log-transform the variables, for use-cases such as having high variance in larger values. Defaults to natural logarithm. See \code{\link[base:Log]{log}} for details and optional args.}
#' \item{'rank'}{rank-transform the variables. Defalts to breaking ties with the average rank of the tied values. See \code{\link[base]{rank}} for details and optional args.}
#' \item{c(opt1, opt2, etc.)}{apply the transform specified in opt1, followed by opt2, etc.}
#' }
#' @param xfm.opts optional arguments to pass to the \code{xfm} option specified. Should be a numbered list of lists, where \code{xfm.opts[[i]]} corresponds to the optional arguments for \code{xfm[i]}. Defaults to the default options for each transform scheme.
#' @param robust whether to use a robust estimate of the covariance matrix when taking PCA. Defaults to \code{FALSE}.
#' @return the svd of X.
#' @author Eric Bridgeford
lol.utils.decomp <- function(X, xfm=FALSE, xfm.opts=list(), ncomp=0, t=.05, robust=FALSE) {
  n <- nrow(X)
  d <- ncol(X)
  # scale if desired before taking SVD
  for (i in 1:length(xfm)) {
    sc <- xfm[i]
    if (!(i %in% names(xfm.opts))) {
      xfm.opts[[i]] <- list()
    }
    if (sc == 'unit') {
      X <- do.call(scale, c(list(X), xfm.opts[[i]]))
    } else if (sc == 'log') {
      X <- do.call(log, c(list(X), xfm.opts[[i]]))
    } else if (sc == 'rank') {
      X <- apply(X, c(2), function(x) {do.call(rank, c(list(x), xfm.opts[[i]]))})
    }
  }
  # take svd
  decomp <- list()
  if (robust) {
    eigenX <- eigen(covRob(X, estim='weighted')$cov)
    decomp$comp <- eigenX$vectors[, 1:ncomp, drop=FALSE]
    decomp$val <- eigenX$values[1:ncomp]
  } else if (ncomp > t*d | ncomp >= d) {
    svdX <- svd(X, nu=0, nv=ncomp)
    decomp$comp <- svdX$v
    decomp$val <- svdX$d
  } else {
    svdX <- irlba(X, nu=0, nv=ncomp)
    decomp$comp <- svdX$v
    decomp$val <- svdX$d
  }
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
#' \item{'unit'}{unit transform the variables, defaulting to centering and scaling to mean 0, variance 1. See \code{\link[base]{scale}} for details and optional args.}
#' \item{'log'}{log-transform the variables, for use-cases such as having high variance in larger values. Defaults to natural logarithm. See \code{\link[base:Log]{log}} for details and optional args.}
#' \item{'rank'}{rank-transform the variables. Defalts to breaking ties with the average rank of the tied values. See \code{\link[base]{rank}} for details and optional args.}
#' \item{c(opt1, opt2, etc.)}{apply the transform specified in opt1, followed by opt2, etc.}
#' }
#' @param xfm.opts optional arguments to pass to the \code{xfm} option specified. Should be a numbered list of lists, where \code{xfm.opts[[i]]} corresponds to the optional arguments for \code{xfm[i]}. Defaults to the default options for each transform scheme.
#' @param robust whether to use a robust estimate of the covariance matrix when taking PCA. Defaults to \code{FALSE}.
#' @param ... trailing args.
#' @return A list containing the following:
#' \item{\code{A}}{\code{[d, r]} the projection matrix from \code{d} to \code{r} dimensions.}
#' \item{\code{d}}{the eigen values associated with the eigendecomposition.}
#' \item{\code{ylabs}}{\code{[K]} vector containing the \code{K} unique, ordered class labels.}
#' \item{\code{centroids}}{\code{[K, d]} centroid matrix of the \code{K} unique, ordered classes in native \code{d} dimensions.}
#' \item{\code{priors}}{\code{[K]} vector containing the \code{K} prior probabilities for the unique, ordered classes.}
#' \item{\code{Xr}}{\code{[n, r]} the \code{n} data points in reduced dimensionality \code{r}.}
#' \item{\code{cr}}{\code{[K, r]} the \code{K} centroids in reduced dimensionality \code{r}.}
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
#' model <- lol.project.lrlda(X=X, Y=Y, r=2)  # use lrlda to project into 2 dimensions
#' @export
lol.project.lrlda <- function(X, Y, r, xfm=FALSE, xfm.opts=list(), robust=FALSE, ...) {
  # class data
  classdat <- lol.utils.info(X, Y, robust=robust)
  priors <- classdat$priors; centroids <- t(classdat$centroids)
  K <- classdat$K; ylabs <- classdat$ylabs
  n <- classdat$n; d <- classdat$d
  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }

  # subtract column means per-class
  Yidx <- sapply(Y, function(y) which(ylabs == y))
  # form class-conditional data matrix
  Xt <- X - centroids[Yidx,]
  # compute the standard projection but with the pre-centered data.
  X.decomp <- lol.utils.decomp(Xt, xfm=xfm, xfm.opts=xfm.opts, ncomp=r, robust=robust)

  return(list(A=X.decomp$comp, d=X.decomp$val, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=lol.embed(X, X.decomp$comp), cr=lol.embed(centroids, X.decomp$comp)))
}
