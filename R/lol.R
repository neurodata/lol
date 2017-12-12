#' Linear Optimal Low-Rank Projection (LOL)
#'
#' A function for implementing the Linear Optimal Low-Rank Projection (LOL) Algorithm.
#'
#' @import irlba
#' @param X [n, d] the data with n samples in d dimensions.
#' @param Y [n] the labels of the samples.
#' @param r the rank of the projection.
#' @return A [d, r] the data projected into r dimensions.
#' @author Eric Bridgeford, adapted from Joshua Vogelstein
#' @export
fs.project.lol <- function(X, Y, r) {
  ylabs <- sort(unique(Y))
  C <- length(ylabs)
  n <- length(Y)
  d <- dim(X)[2]
  nv <- r - C

  A <- array(0, dim=c(d, r))
  for (i in 1:length(ylabs)) {
    select <- Y == ylabs[i]
    A[,i] <- colMeans(X[select,,drop=FALSE])
  }

  svd <- irlba(t(as.matrix(X)), nv=0, nu=nv)
  A[,C:r] <- svd$u

  # orthogonalize and normalize
  A <- qr.Q(qr(A))
  return(A)
}
