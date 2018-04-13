#' Bayes Optimal
#'
#' A function for recovering the Bayes Optimal Projection, which optimizes Bayes classification.
#'
#' @import irlba
#' @importFrom MASS ginv
#' @param X \code{[n, p]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param mus \code{[d, K]} the \code{K} class means in \code{d} dimensions.
#' @param Sigmas \code{[d, d, K]} the \code{K} class covariance matrices in \code{d} dimensions.
#' @param priors \code{[K]} the priors for each of the \code{K} classes.
#' @param ... optional args.
#' @return A list of class \code{embedding} containing the following:
#' \item{\code{A}}{\code{[d, r]} the projection matrix from \code{d} to \code{r} dimensions.}
#' \item{\code{d}}{the eigen values associated with the eigendecomposition.}
#' \item{\code{ylabs}}{\code{[K]} vector containing the \code{K} unique, ordered class labels.}
#' \item{\code{centroids}}{\code{[K, d]} centroid matrix of the \code{K} unique, ordered classes in native \code{d} dimensions.}
#' \item{\code{priors}}{\code{[K]} vector containing the \code{K} prior probabilities for the unique, ordered classes.}
#' \item{\code{Xr}}{\code{[n, r]} the \code{n} data points in reduced dimensionality \code{r}.}
#' \item{\code{cr}}{\code{[K, r]} the \code{K} centroids in reduced dimensionality \code{r}.}
#' @author Eric Bridgeford
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' # obtain bayes-optimal projection of the data
#' model <- lol.project.bayes_optimal(X=X, Y=Y, mus=data$mus,
#'                                    S=data$Sigmas, priors=data$priors)
#' @export
lol.project.bayes_optimal <- function(X, Y, mus, Sigmas, priors, ...) {
  info <- lol.utils.info(X, Y)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  deltas <- lol.utils.deltas(centroids, priors)
  centroids <- t(centroids)
  E <- lol.mvr(Sigmas, mus, priors)

  A <- ginv(E) %*% deltas
  return(structure(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
                        Xr=lol.embed(X, A), cr=lol.embed(centroids, A)), class="embedding"))
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
