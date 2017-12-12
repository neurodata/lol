#' Low-Rank Linear Discriminant Analysis (lrlda)
#'
#' A function for implementing the Low-Rank Linear Discriminant Analysis (lrlda) Projection.
#'
#' @import irlba
#' @import MASS
#' @param X [n, d] the data with n samples in d dimensions.
#' @param Y [n] the labels of the samples, with C classes.
#' @return W [d, r] the projection matrix for the linearly optimal projection from d to C-1 dimensions.
#' @author Ran Liu
#' @export
fs.project.lrlda <- function(X, Y) {
  ylabs <- sort(unique(Y))
  C <- length(ylabs)
  n <- length(Y)
  d <- dim(X)[2]
  # Compute prior probabilities
  pi <- sapply(ylabs, function(y) sum(Y==y)/length(Y))
  counts <- sapply(ylabs, function(y) sum(Y==y))
  # Compute class means
  class_means <- sapply(ylabs, function(y) colMeans(X[Y==y,,drop=FALSE]))
  overall_mean <- colMeans(X)

  sb <- array(0, dim=c(d, d))
  sw <- array(0, dim=c(d, d))

  for (i in 1:length(ylabs)) {
    sb <- sb + counts[i] * ((class_means[i,] - overall_mean) %*% t(class_means[i,] - overall_mean))
    x <- X[Y==ylabs[i],]
    sw <- sw + t(x-class_means[i,]) %*% (x-class_means[i,])
  }

  A <- eigen(ginv(sw) %*% sb)
  return(A$vectors[,1:length(ylabs)])
}
