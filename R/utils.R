#' A function that performs basic utilities about the data.
#' @importFrom robustbase colMedians
#' @param X \code{[n, d]} the data with n samples in d dimensions.
#' @param Y \code{[n]} the labels of the samples.
#' @param robust whether to perform PCA on a robust estimate of the covariance matrix or not. Defaults to \code{FALSE}.
#' @param ... optional args.
#' @return \code{n} the number of samples.
#' @return \code{d} the number of dimensions.
#' @return ylabs \code{[K]} vector containing the unique, ordered class labels.
#' @return priors \code{[K]} vector containing prior probability for the unique, ordered classes.
#' @author Eric Bridgeford
lol.utils.info <- function(X, Y, robust=FALSE, ...) {
  ylabs <- as.vector(sort(unique(Y)))
  dimx <- dim(X)
  n <- dimx[1]; d <- dimx[2]
  K <- length(ylabs)
  # compute the fraction of each class for prior p
  priors = sapply(ylabs, function(y) sum(Y == y)/n)
  if (!robust) {
    centroids <- as.matrix(array(sapply(ylabs, function(y) colMeans(X[Y==y,,drop=FALSE])), dim=c(d, K)))
  } else {
    # robust estimator of the mean is the median
    centroids <- as.matrix(array(sapply(ylabs, function(y) colMedians(X[Y==y,,drop=FALSE])), dim=c(d, K)))
  }
  return(list(n=n, d=d, ylabs=ylabs, priors=priors, centroids=centroids, K=K))
}

#' A function that performs a utility computation of information about the differences of the classes.
#' @param centroids \code{[d, K]} centroid matrix of the unique, ordered classes.
#' @param priors \code{[K]} vector containing prior probability for the unique, ordered classes.
#' @param ... optional args.
#' @return deltas \code{[d, K]} the K difference vectors.
#' @author Eric Bridgeford
lol.utils.deltas <- function(centroids, priors, ...) {
  d <- dim(centroids)[1]; K <- length(priors)
  # compute the rank-K difference space as deltas(i) = mus(i) - mus(0) where the mus are ordered
  # by decreasing prior
  deltas <- array(0, dim=c(d, K))
  srt_prior <- sort(priors, decreasing=TRUE, index.return=TRUE)$ix
  deltas[,1] <- centroids[,srt_prior[1]]
  for (i in 2:K) {
    deltas[,i] <- deltas[,1] - centroids[,srt_prior[2]]
  }
  return(deltas[,2:K,drop=FALSE])
}

#' A function for one-hot encoding categorical respose vectors.
#' @param Y [n] a vector of the categorical resposes, with \code{K} unique categories.
#' @return a list containing the following:
#' \item{Yh}{[n, K] the one-hot encoded Y respose variable.}
#' \item{ylabs}{[K] a vector of the y names corresponding to each response column.}
#' @author Eric Bridgeford
lol.utils.ohe <- function(Y) {
  ylabs <- unique(Y)
  K <- length(ylabs)
  n <- length(Y)
  # hot-encoding of Y categorical variables
  Yh <- array(0, dim=c(n, K))
  # Yind is a indicator of membership in each respective class
  for (i in 1:K) {
    Yh[Y == ylabs[i],i] <- 1
  }
  colnames(Yh) = ylabs
  return(Yh)
}
