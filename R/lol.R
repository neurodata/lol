#' Linear Optimal Low-Rank Projection (LOL)
#'
#' A function for implementing the Linear Optimal Low-Rank Projection (LOL) Algorithm.
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param r the rank of the projection. Note that \code{r >= K}, and \code{r < d}.
#' @param xfm whether to transform the variables before taking the SVD.
#' \itemize{
#' \item{\code{FALSE} apply no transform to the variables.}
#' \item{\code{'unit'} unit transform the variables, defaulting to centering and scaling to mean 0, variance 1. See \link[base]{scale} for details and optional args.}
#' \item{\code{'log'} log-transform the variables, for use-cases such as having high variance in larger values. Defaults to natural logarithm. See \link[base]{log} for details and optional args.}
#' \item{\code{'rank'} rank-transform the variables. Defalts to breaking ties with the average rank of the tied values. See \link[base]{rank} for details and optional args.}
#' \item{\code{c(opt1, opt2, etc.)} apply the transform specified in opt1, followed by opt2, etc.}
#' }
#' @param xfm.opts optional arguments to pass to the \code{xfm} option specified. Should be a numbered list of lists, where \code{xfm.opts[[i]]} corresponds to the optional arguments for \code{xfm[i]}. Defaults to the default options for each transform scheme.
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
#' @export
lol.project.lol <- function(X, Y, r, xfm=FALSE, xfm.opts=list(), ...) {
  # class data
  info <- lol.utils.info(X, Y)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }
  deltas <- lol.utils.deltas(centroids, priors)
  centroids <- t(centroids)

  nv <- r - (K)
  lrlda <- list(d=NULL)
  if (nv > 0) {
    lrlda <- lol.project.lrlda(X, Y, r=nv, xfm=xfm, xfm.opts=xfm.opts)
    A <- cbind(deltas, lrlda$A)
  } else {
    A <- deltas[, 1:r, drop=FALSE]
  }

  # orthogonalize and normalize
  A <- qr.Q(qr(A))
  return(list(A=A, d=lrlda$d, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=lol.embed(X, A), cr=lol.embed(centroids, A)))
}

#' Quadratic Optimal QDA (QOQ)
#'
#' A function for implementing the Quadratic Optimal QDA Projection (QOQ) Algorithm, an intuitive adaptation of the Linear Optimal Low-Rank Projection (LOL).
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param r the rank of the projection. Note that \code{r >= K}, and \code{r < d}.
#' @param xfm whether to transform the variables before taking the SVD.
#' \itemize{
#' \item \code{FALSE} apply no transform to the variables.
#' \item \code{'unit'} unit transform the variables, defaulting to centering and scaling to mean 0, variance 1. See \link[base]{scale} for details and optional args.
#' \item \code{'log'} log-transform the variables, for use-cases such as having high variance in larger values. Defaults to natural logarithm. See \link[base]{log} for details and optional args.
#' \item \code{'rank'} rank-transform the variables. Defalts to breaking ties with the average rank of the tied values. See \link[base]{rank} for details and optional args.
#' \item \code{c(opt1, opt2, etc.)} apply the transform specified in opt1, followed by opt2, etc.
#' }
#' @param xfm.opts optional arguments to pass to the \code{xfm} option specified. Should be a numbered list of lists, where \code{xfm.opts[[i]]} corresponds to the optional arguments for \code{xfm[i]}. Defaults to the default options for each transform scheme.
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
#' \code{vignette("qoq", package = "lolR")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(lolR)
#' data <- lol.sims.qdtoep(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
  #' model <- lol.project.qoq(X=X, Y=Y, r=5)  # use qoq to project into 5 dimensions
#' @export
lol.project.qoq <- function(X, Y, r, xfm=FALSE, xfm.opts=list(), ...) {
  # class data
  info <- lol.utils.info(X, Y)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }
  deltas <- lol.utils.deltas(centroids, priors)
  centroids <- t(centroids)

  nv <- r - (K)
  Aclass <- array(0, dim=c(d, 0))  # the class-wise egvecs
  vclass <- c()  # the class-wise egvals
  vclass.res <- list(d=NULL)
  if (nv > 0) {
    for (ylab in ylabs) {
      Xclass = X[Y == ylab,]
      obj <- lol.project.pca(Xclass, r=nv, xfm=xfm, xfm.opts=xfm.opts)
      Aclass <- cbind(Aclass, obj$A)
      vclass <- c(vclass, obj$d[1:nv])
    }
    # take the nv from the A computed for each class using the
    # nv with the top eigenvalues from Aclass
    A <- cbind(deltas, Aclass[, sort(vclass, index.return=TRUE)$ix[1:nv]])
    vclass.res$d <- sort(vclass)[1:nv]
  } else {
    A <- deltas[, 1:r, drop=FALSE]
  }

  # orthogonalize and normalize
  A <- qr.Q(qr(A))
  return(list(A=A, d=vclass.res$d, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=lol.embed(X, A), cr=lol.embed(centroids, A)))
}

#' Partial Least Squares Optimal Low-Rank Projection (PLSOL)
#'
#' A function for implementing the Partial Least Squares Optimal Low-Rank Projection Projection (PLSOL) Algorithm.
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param r the rank of the projection. Note that \code{r >= K}, and \code{r < d}.
#' @param xfm whether to transform the variables before taking the SVD.
#' \itemize{
#' \item{\code{FALSE} apply no transform to the variables.}
#' \item{\code{'unit'} unit transform the variables, defaulting to centering and scaling to mean 0, variance 1. See \link[base]{scale} for details and optional args.}
#' \item{\code{'log'} log-transform the variables, for use-cases such as having high variance in larger values. Defaults to natural logarithm. See \link[base]{log} for details and optional args.}
#' \item{\code{'rank'} rank-transform the variables. Defalts to breaking ties with the average rank of the tied values. See \link[base]{rank} for details and optional args.}
#' \item{\code{c(opt1, opt2, etc.)} apply the transform specified in opt1, followed by opt2, etc.}
#' }
#' @param xfm.opts optional arguments to pass to the \code{xfm} option specified. Should be a numbered list of lists, where \code{xfm.opts[[i]]} corresponds to the optional arguments for \code{xfm[i]}. Defaults to the default options for each transform scheme.
#' @param ... trailing args.
#' @return A list containing the following:
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
#' model <- lol.project.plsol(X=X, Y=Y, r=5)  # use lol to project into 5 dimensions
#' @export
lol.project.plsol <- function(X, Y, r, xfm=FALSE, xfm.opts=list(), ...) {
  # class data
  info <- lol.utils.info(X, Y)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }
  centroids <- t(centroids)


  A.pls <- lol.project.pls(X, Y, r=min(r, K-1))$A

  nv <- r - min(r, (K - 1))
  lrlda <- list(d=NULL)
  if (nv > 0) {
    lrlda <- lol.project.lrlda(X, Y, r=nv, xfm=xfm, xfm.opts=xfm.opts)
    A <- cbind(A.pls, lrlda$A)
  } else {
    A <- A.pls[, 1:r, drop=FALSE]
  }

  # orthogonalize and normalize
  A <- qr.Q(qr(A))
  return(list(A=A, d=lrlda$d, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=lol.embed(X, A), cr=lol.embed(centroids, A)))
}
#' Partial Least Squares Optimal Low-Rank Projection K (PLSOLK)
#'
#' A function for implementing the Partial Least Squares Optimal Low-Rank Projection Projection K (PLSOLK) Algorithm.
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param r the rank of the projection. Note that \code{r >= K}, and \code{r < d}.
#' @param xfm whether to transform the variables before taking the SVD.
#' \itemize{
#' \item{\code{FALSE} apply no transform to the variables.}
#' \item{\code{'unit'} unit transform the variables, defaulting to centering and scaling to mean 0, variance 1. See \link[base]{scale} for details and optional args.}
#' \item{\code{'log'} log-transform the variables, for use-cases such as having high variance in larger values. Defaults to natural logarithm. See \link[base]{log} for details and optional args.}
#' \item{\code{'rank'} rank-transform the variables. Defalts to breaking ties with the average rank of the tied values. See \link[base]{rank} for details and optional args.}
#' \item{\code{c(opt1, opt2, etc.)} apply the transform specified in opt1, followed by opt2, etc.}
#' }
#' @param xfm.opts optional arguments to pass to the \code{xfm} option specified. Should be a numbered list of lists, where \code{xfm.opts[[i]]} corresponds to the optional arguments for \code{xfm[i]}. Defaults to the default options for each transform scheme.
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
#'
#' @author Eric Bridgeford
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.plsolk(X=X, Y=Y, r=5)  # use lol to project into 5 dimensions
#' @export
lol.project.plsolk <- function(X, Y, r, xfm=FALSE, xfm.opts=list(), ...) {
  # class data
  info <- lol.utils.info(X, Y)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }
  A.pls <- lol.project.pls(X, Y, r=min(r, K))$A
  centroids <- t(centroids)

  nv <- r - min(r, (K))
  lrlda <- list(d=NULL)
  if (nv > 0) {
    lrlda <- lol.project.lrlda(X, Y, r=nv, xfm=xfm, xfm.opts=xfm.opts)
    A <- cbind(A.pls, lrlda$A)
  } else {
    A <- A.pls[, 1:r, drop=FALSE]
  }

  # orthogonalize and normalize
  A <- qr.Q(qr(A))
  return(list(A=A, d=lrlda$d, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=lol.embed(X, A), cr=lol.embed(centroids, A)))
}
