#' A function that performs a utility computation of prior information about the data.
#' @import irlba
#' @param X [n, d] the data with n samples in d dimensions.
#' @param Y [n] the labels of the samples.
#' @return ylabs [C] vector containing the unique, ordered class labels.
#' @return centroids [C, d] centroid matrix of the unique, ordered classes.
#' @return priors [C] vector containing prior probability for the unique, ordered classes.
#' @return n the number of samples.
#' @return d the number of dimensions.
#' @author Eric Bridgeford
gs.utils.classdat <- function(X, Y) {
  ylabs <- sort(unique(Y))
  dimx <- dim(X)
  n <- dimx[1]; d <- dimx[2]

  # compute the fraction of each class for prior p
  priors = sapply(ylabs, function(y) sum(Y == y)/n)
  # compute the centroids of each class
  centroids <- t(as.matrix(sapply(ylabs, function(y) colMeans(X[Y==y,,drop=FALSE]))))

  return(list(ylabs=ylabs, centroids=centroids, priors=priors, C=length(ylabs),
              n=n, d=d))
}
