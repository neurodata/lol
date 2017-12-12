#' Linear Optimal Low-Rank Projection (LOL)
#'
#' A function for implementing the Linear Optimal Low-Rank Projection (LOL) Algorithm.
#'
#' @import irlba
#' @param X [n, d] the data with n samples in d dimensions.
#' @param Y [n] the labels of the samples.
#' @param r the rank of the projection.
#' @return A [d, r] the projection matrix for the linearly optimal projection from d to r dimensions.
#' @author Eric Bridgeford
#' @export
fs.project.lol <- function(X, Y, r) {
  ylabs <- sort(unique(Y))
  C <- length(ylabs)
  n <- length(Y)
  d <- dim(X)[2]
  nv <- r - C

  A <- sapply(ylabs, function(y) colMeans(X[Y==y,,drop=FALSE]))

  svd <- irlba(t(as.matrix(X)), nv=0, nu=nv)
  A <- cbind(A, svd$u)

  # orthogonalize and normalize
  A <- qr.Q(qr(A))
  return(A)
}
