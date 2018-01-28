#' Linear Optimal Low-Rank Projection (LOL)
#'
#' A function for implementing the Linear Optimal Low-Rank Projection (LOL) Algorithm.
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param r the rank of the projection. Note that \code{r >= K}, and \code{r < d}.
#' @param ... optional args.
#' @return A list of class \code{embedding} containing the following:
#' \item{A}{\code{[d, r]} the projection matrix from \code{d} to \code{r} dimensions.}
#' \item{ylabs}{\code{[K]} vector containing the \code{K} unique, ordered class labels.}
#' \item{centroids}{\code{[K, d]} centroid matrix of the \code{K} unique, ordered classes in native \code{d} dimensions.}
#' \item{priors}{\code{[K]} vector containing the \code{K} prior probabilities for the unique, ordered classes.}
#' \item{Xr}{\code{[n, r]} the \code{n} data points in reduced dimensionality \code{r}.}
#' \item{cr}{\code{[K, r]} the \code{K} centroids in reduced dimensionality \code{r}.}
#' @author Eric Bridgeford
#' @examples
#' library(lol)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.lol(X=X, Y=Y, r=5)  # use lol to project into 5 dimensions
#' @export
lol.project.lol <- function(X, Y, r, ...) {
  # class data
  info <- lol.utils.info(X, Y)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  deltas <- lol.utils.deltas(centroids, priors)
  centroids <- t(centroids)

  nv <- r - (K)
  if (nv > 0) {
    A <- cbind(deltas, lol.project.cpca(X, Y, nv)$A)
  } else {
    A <- deltas
    A <- A[,1:r,drop=FALSE]
  }

  # orthogonalize and normalize
  A <- qr.Q(qr(A))
  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=lol.embed(X, A), cr=lol.embed(centroids, A)))
}
