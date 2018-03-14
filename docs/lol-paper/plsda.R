require(caret)
require(pls)
require(R.utils)

lol.project.pls <- function(X, Y, r, ...) {
  info <- lol:::lol.utils.info(X, Y)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  A <- plsda(X, as.factor(Y), r)$projection
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
    plsA <- plsda(X, as.factor(Y), r)$projection
    A <- cbind(deltas, plsA)
  } else {
    A <- deltas[, 1:r, drop=FALSE]
  }
  # orthogonalize and normalize
  A <- qr.Q(qr(A))
  return(list(A=A, Xa=lol.embed(X, A), ylabs=ylabs))
}

lol.project.opals <- function(X, Y, r, ...) {
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
    plsA <- plsda(Xt, as.factor(Y), r)$projection
    A <- cbind(deltas, plsA)
  } else {
    A <- deltas[, 1:r, drop=FALSE]
  }
  # orthogonalize and normalize
  A <- qr.Q(qr(A))
  return(list(A=A, Xa=lol.embed(X, A), ylabs=ylabs))
}
