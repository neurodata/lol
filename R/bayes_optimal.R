#' Bayes Optimal
#'
#' A function for recovering the Bayes Optimal Projection, which optimizes Bayes classification.
#'
#' @import irlba
#' @param X [n, p] the data with n samples in d dimensions.
#' @param mus [d, K] the true means of each class.
#' @param Sigmas [d, d, K] the true covariances of each class.
#' @param priors [K] the priors for each class.
#' @return A [d, r] the projection matrix from d to r dimensions.
#' @return ylabs [K] vector containing the unique, ordered class labels.
#' @return centroids [K, d] centroid matrix of the unique, ordered classes.
#' @return priors [K] vector containing prior probability for the unique, ordered classes.
#' @return Xr [n, r] the data in reduced dimensionality.
#' @return cr [K, r] the centroids in reduced dimensionality.
#' @author Jason Yim and Eric Bridgeford
#' @examples
#' library(lol)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.bayes_optimal(X=X, mus=data$mus, S=data$Sigmas, priors=data$priors)  # obtain bayes-optimal projection of the data
#' @export
lol.project.bayes_optimal <- function(X, Y, mus, Sigmas, priors, ...) {
  info <- lol:::lol.utils.info(X, Y)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  deltas <- lol:::lol.utils.deltas(centroids, priors)
  centroids <- t(centroids)
  E <- lol.mvr(Sigmas, mus, priors)

  A <- MASS::ginv(E) %*% deltas
  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=X %*% A, cr=centroids %*% A))
}

# a function to compute the joint-variance analytically given
# data simulated from multiple gaussians
lol.mvr <- function(S, mus, priors) {
  K <- dim(mus)[2]
  d <- dim(mus)[1]
  muh <- apply(abind(lapply(1:K, function(i) {
    mus[,i]*priors[i]
  }), along=2), c(1), sum)
  T1 <- apply(abind(lapply(1:K, function(i) {
    priors[i]*S[,,i]
  }), along=3), c(1,2), sum)
  T2 <- apply(abind(lapply(1:K, function(i) {
    priors[i]*(mus[,i,drop=FALSE] - muh) %*% t(mus[,i,drop=FALSE])
  }), along=3), c(1,2), sum)
  Sa <- T1 + T2
  return(Sa)
}
