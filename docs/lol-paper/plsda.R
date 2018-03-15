require(pls)
require(R.utils)

lol.project.pls <- function(X, Y, r, ...) {
  info <- lol:::lol.utils.info(X, Y)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d

  Yh <- lol:::lol.utils.ohe(Y)
  tmpData <- data.frame(n=paste("row", 1:nrow(Yh), sep=""))
  tmpData$Y <- Yh
  tmpData$X <- X

  A <- pls::plsr(Y ~ X, data=tmpData, ncomp=r)$projection
  return(list(A=A, Xa=lol.embed(X, A), ylabs=ylabs))
}

lol.project.mpls <- function(X, Y, r, ...) {
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
  if (nv > 0) {
    # compute the standard projection but with the pre-centered data.
    plsA <- lol.project.pls(X, Y, r)$A
    A <- cbind(deltas, plsA)
  } else {
    A <- deltas[, 1:r, drop=FALSE]
  }
  # orthogonalize and normalize
  A <- qr.Q(qr(A))
  return(list(A=A, Xa=lol.embed(X, A), ylabs=ylabs))
}

lol.project.opal <- function(X, Y, r, ...) {
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
  if (nv > 0) {
    # subtract column means per-class
    Yidx <- sapply(Y, function(y) which(ylabs == y))
    # form class-conditional data matrix
    Xt <- X - centroids[Yidx,]
    # compute the standard projection but with the pre-centered data.
    plsA <- lol.project.pls(Xt, Y, r)$A
    A <- cbind(deltas, plsA)
  } else {
    A <- deltas[, 1:r, drop=FALSE]
  }
  # orthogonalize and normalize
  A <- qr.Q(qr(A))
  return(list(A=A, Xa=lol.embed(X, A), ylabs=ylabs))
}
