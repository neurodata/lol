#' Low-rank Canonical Correlation Analysis (LR-CCA)
#'
#' A function for implementing the Low-rank Canonical Correlation Analysis (LR-CCA) Algorithm.
#'
#' @import irlba
#' @param X [n, p] the data with n samples in d dimensions.
#' @param Y [n, q] the labels of the samples.
#' @param r the rank of the projection.
#' @return A [d, r] the projection matrix for the linearly optimal projection from d to r dimensions.
#' @author jason Yim
#' @export
low_rank_cca <- function(X, Y, r) {
  cxy <- cancor(X,Y)
  return(cxy$xcoef[,rank(cxy$cor)[1:r]][,1:r,drop=FALSE])
}
