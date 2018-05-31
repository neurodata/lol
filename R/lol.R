#' Linear Optimal Low-Rank Projection (LOL)
#'
#' A function for implementing the Linear Optimal Low-Rank Projection (LOL) Algorithm. This algorithm allows users to find an optimal
#' projection from `d` to `r` dimensions, where `r << d`, by combining information from the first and second moments in thet data.
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param r the rank of the projection. Note that \code{r >= K}, and \code{r < d}.
#' @param first.moment the function to capture the first moment. Defaults to \code{'delta'}.
#' \itemize{
#' \item{\code{'delta'} capture the first moment with the hyperplane separating the per-class means.}
#' \item{\code{FALSE} do not capture the first moment.}
#' }
#' @param second.moment the function to capture the second moment. Defaults to \code{'linear'}.
#' \itemize{
#' \item{\code{'linear'} performs PCA on the class-conditional data to capture the second moment, retaining the vectors with the top singular values.}
#' \item{\code{'quadratic'} performs PCA on the data for each class separately to capture the second moment, retaining the vectors with the top singular values from each class's PCA.}
#' \item{\code{'pls'} performs PLS on the data to capture the second moment, retaining the vectors that maximize the correlation between the different classes.}
#' \item{\code{FALSE} do not capture the second moment.}
#' }
#' @param orthogonalize whether to orthogonalize the projection matrix. Defaults to \code{FALSE}.
#' @param second.moment.xfm whether to transform the variables before taking the SVD.
#' \itemize{
#' \item{\code{FALSE} apply no transform to the variables.}
#' \item{\code{'unit'} unit transform the variables, defaulting to centering and scaling to mean 0, variance 1. See \link[base]{scale} for details and optional args.}
#' \item{\code{'log'} log-transform the variables, for use-cases such as having high variance in larger values. Defaults to natural logarithm. See \link[base]{log} for details and optional args.}
#' \item{\code{'rank'} rank-transform the variables. Defalts to breaking ties with the average rank of the tied values. See \link[base]{rank} for details and optional args.}
#' \item{\code{c(opt1, opt2, etc.)} apply the transform specified in opt1, followed by opt2, etc.}
#' }
#' @param second.moment.xfm.opts optional arguments to pass to the \code{xfm} option specified. Should be a numbered list of lists, where \code{second.moment.xfm.opts[[i]]} corresponds to the optional arguments for \code{second.moment.xfm[[i]]}.
#' Defaults to the default options for each transform scheme.
#' @param ... trailing args.
#' @return A list containing the following:
#' \item{\code{A}}{\code{[d, r]} the projection matrix from \code{d} to \code{r} dimensions.}
#' \item{\code{d}}{the eigen values associated with the eigendecomposition.}
#' \item{\code{ylabs}}{\code{[K]} vector containing the \code{K} unique, ordered class labels.}
#' \item{\code{centroids}}{\code{[K, d]} centroid matrix of the \code{K} unique, ordered classes in native \code{d} dimensions.}
#' \item{\code{priors}}{\code{[K]} vector containing the \code{K} prior probabilities for the unique, ordered classes.}
#' \item{\code{Xr}}{\code{[n, r]} the \code{n} data points in reduced dimensionality \code{r}.}
#' \item{\code{cr}}{\code{[K, r]} the \code{K} centroids in reduced dimensionality \code{r}.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("lol", package = "lolR")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.lol(X=X, Y=Y, r=5)  # use lol to project into 5 dimensions
#'
#' model <- lol.project.lol(X=X, Y=Y, r=5, second.moment='quadratic')  # use qoq to project into 5 dimensions
#' @export
lol.project.lol <- function(X, Y, r, second.moment.xfm=FALSE, second.moment.xfm.opts=list(),
                            first.moment='delta', second.moment='linear', orthogonalize=FALSE, ...) {
  # class data
  info <- lol.utils.info(X, Y)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }

  if (first.moment == "delta") {
    first.moment.proj <- lol.utils.deltas(centroids, priors)[, 2:K, drop=FALSE]
  } else {
    first.moment.proj <- array(0, dim=c(d, 0))
  }

  nv <- r - dim(first.moment.proj)[2]

  if (second.moment == "linear" & nv > 0) {
    lrlda <- lol.project.lrlda(X, Y, r=nv, xfm=second.moment.xfm, xfm.opts=second.moment.xfm.opts)
    #d <- lrlda$d
    second.moment.proj <- lrlda$A
  } else if (second.moment == "quadratic" & nv > 0) {
    Aclass <- array(0, dim=c(d, 0))  # the class-wise egvecs
    vclass <- c()  # the class-wise egvals
    for (ylab in ylabs) {
      Xclass = X[Y == ylab,]
      obj <- lol.project.pca(Xclass, r=nv, xfm=second.moment.xfm, xfm.opts=second.moment.xfm.opts)
      Aclass <- cbind(Aclass, obj$A)
      vclass <- c(vclass, obj$d[1:nv])
    }
    # take the nv from the A computed for each class using the
    # nv with the top eigenvalues from Aclass
    second.moment.proj <- Aclass[, sort(vclass, index.return=TRUE, decreasing=TRUE)$ix[1:nv]]
    #d <- sort(vclass)[1:nv]
  } else if (second.moment == "pls" & nv > 0) {
    pls.res <- lol.project.pls(X, Y, nv)
    second.moment.proj <- pls.res$A
  } else {
    second.moment.proj <- array(0, dim=c(d, 0))
  }

  # combine estimates of first and second moment
  A <- cbind(first.moment.proj, second.moment.proj)
  if (orthogonalize) {
    A <- qr.Q(qr(A))
  }
  centroids <- t(centroids)
  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=lol.embed(X, A), cr=lol.embed(centroids, A), second.moment=second.moment,
              first.moment=first.moment))
}
