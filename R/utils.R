#' A function that performs basic utilities about the data.
#' @import irlba
#' @param X \code{[n, d]} the data with n samples in d dimensions.
#' @param Y \code{[n]} the labels of the samples.
#' @param ... optional args.
#' @return \code{n} the number of samples.
#' @return \code{d} the number of dimensions.
#' @return ylabs \code{[K]} vector containing the unique, ordered class labels.
#' @return priors \code{[K]} vector containing prior probability for the unique, ordered classes.
#' @author Eric Bridgeford
lol.utils.info <- function(X, Y, ...) {
  ylabs <- as.vector(sort(unique(Y)))
  dimx <- dim(X)
  n <- dimx[1]; d <- dimx[2]
  K <- length(ylabs)
  # compute the fraction of each class for prior p
  priors = sapply(ylabs, function(y) sum(Y == y)/n)
  centroids <- as.matrix(sapply(ylabs, function(y) colMeans(X[Y==y,,drop=FALSE])))
  return(list(n=n, d=d, ylabs=ylabs, priors=priors, centroids=centroids, K=K))
}

#' A function that performs a utility computation of information about the differences of the classes.
#' @import irlba
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
  gr_mix <- srt_prior[1]
  deltas[,1] <- centroids[,gr_mix]
  for (i in 2:K) {
    deltas[,(i)] <- centroids[,srt_prior[i]] - deltas[,1]
  }
  return(deltas)
}
