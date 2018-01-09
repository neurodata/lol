#' Linear Optimal Low-Rank Projection (LOL)
#'
#' A function for implementing the Linear Optimal Low-Rank Projection (LOL) Algorithm.
#'
#' @param X [n, d] the data with n samples in d dimensions.
#' @param Y [n] the labels of the samples.
#' @param r the rank of the projection. Note that r >= length(unique(Y)).
#' @return A [d, r] the projection matrix from d to r dimensions.
#' @return ylabs [K] vector containing the unique, ordered class labels.
#' @return centroids [K, d] centroid matrix of the unique, ordered classes.
#' @return priors [K] vector containing prior probability for the unique, ordered classes.
#' @return Xr [n, r] the data in reduced dimensionality.
#' @return cr [K, r] the centroids in reduced dimensionality.
#' @author Eric Bridgeford
#' @examples
#' library(fselect)
#' data <- fs.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- fs.project.lol(X=X, Y=Y, r=5)  # use lol to project into 5 dimensions
#' @export
fs.project.lol <- function(X, Y, r) {
  # class data
  classdat <- fselect:::fs.utils.classdat(X, Y)
  priors <- classdat$priors; centroids <- classdat$centroids
  K <- classdat$K; ylabs <- classdat$ylabs
  n <- classdat$n; d <- classdat$d
  nv <- r - K

  if (r < K) {
    stop(paste("LOL requires the rank of the projection to be >= the number of unique classes in Y.",
               sprintf("rank is r=%d, number of classes is K=%d.", r, K)))
  }

  if (nv > 0) {
    A <- cbind(t(centroids), fs.project.cpca(X, Y, nv)$A)
  } else {
    A <- t(centroids)
  }

  # orthogonalize and normalize
  A <- Matrix::qr.Q(Matrix::qr(A))
  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=X %*% A, cr=centroids %*% A))
}
