require(lolR)
require(FlashR)
require(abind)
require(dplyr)

flashx.pca <- function(X, r, ...) {
  X <- fm.as.matrix(X); Y <- as.vector(Y)
  # mean center by the column mean
  d <- dim(X)[2]
  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }
  # center the data
  Xc  <- sweep(X, 2, colMeans(X), '-')
  X.decomp <- flashx.decomp(Xc, ncomp=r)

  return(list(Xr = flashx.embed(X, X.decomp$comp), d=X.decomp$val, A=X.decomp$comp))
}

flashx.embed <- function(X, A) {
  return(as.matrix(X %*% A))
}

flashx.decomp <- function(X, ncomp=0) {
  svdX <- fm.svd(X, nu=0, nv=ncomp)
  decomp=list(comp=svdX$v, val=svdX$d)
  return(decomp)
}

flashx.lrlda <- function(X, Y, r, ...) {
  X <- fm.as.matrix(X); Y <- as.vector(Y)
  # class data
  if (min(Y) > 0) {
    Y <- Y - min(Y)
  }
  n <- length(Y); d <- ncol(X)
  counts <- as.data.frame(table(fm.conv.FM2R(Y)))
  K <- length(counts$Freq); ylabs <- unique(counts$Var1)
  priors <- counts$Freq/n

  centroids <- fm.mapply.col(
    # compute the group-wise sums
    fm.groupby(X, 2, fm.as.factor(fm.as.vector(Y), K), fm.bo.add),
    # normalize each column by the class mean for that column
    counts$Freq, fm.bo.div)

  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }

  # subtract column means per-class
  Yidx <- sapply(Y, function(y) which(ylabs == y))
  # form class-conditional data matrix
  Xc <- X - fm.get.rows(centroids, Yidx)
  # compute the standard projection but with the pre-centered data
  X.decomp <- flashx.decomp(Xc, ncomp=r)

  return(list(Xr = flashx.embed(X, X.decomp$comp), d=X.decomp$val, A=X.decomp$comp))
}

flashx.deltas <- function(centroids, priors) {
  d <- nrow(centroids); K <- length(priors)
  # compute the rank-K difference space as deltas(i) = mus(i) - mus(0) where the mus are ordered
  # by decreasing prior
  deltas <- list()
  str_prior <- sort(priors, decreasing=TRUE, index.return=TRUE)$ix
  for (i in 2:K) {
    deltas[[i-1]] <- fm.get.rows(centroids, str_prior[1]) - fm.get.rows(centroids, str_prior[i])
  }
  if (length(deltas) > 1) {
    deltas <- fm.cbind.list(deltas)
  } else {
    deltas <- deltas[[1]]
  }
  return(t(deltas))
}

flashx.mdiff <- function(X, Y, ...) {
  X <- fm.as.matrix(X); Y <- as.vector(Y)
  # class data
  if (min(Y) > 0) {
    Y <- Y - min(Y)
  }
  n <- length(Y); d <- ncol(X)
  counts <- as.data.frame(table(fm.conv.FM2R(Y)))
  K <- length(counts$Freq); ylabs <- unique(counts$Var1)
  priors <- counts$Freq/n

  centroids <- fm.mapply.col(
    # compute the group-wise sums
    fm.groupby(X, 2, fm.as.factor(fm.as.vector(Y), K), fm.bo.add),
    # normalize each column by the class mean for that column
    classdat, fm.bo.div)

  deltas <- flashx.deltas(centroids, priors)
  return(list(Xr=flashx.embed(X, deltas), A=deltas))
}

flashx.lol <- function(X, Y, r, ...) {
  X <- fm.as.matrix(X); Y <- as.vector(Y)
  # class data
  if (min(Y) > 0) {
    Y <- Y - min(Y)
  }
  n <- length(Y); d <- ncol(X)
  counts <- as.data.frame(table(fm.conv.FM2R(Y)))
  K <- length(counts$Freq); ylabs <- unique(counts$Var1)
  priors <- counts$Freq/n

  centroids <- fm.mapply.col(
    # compute the group-wise sums
    fm.groupby(X, 2, fm.as.factor(fm.as.vector(Y), K), fm.bo.add),
    # normalize each column by the class mean for that column
    counts$Freq, fm.bo.div)

  deltas <- flashx.deltas(centroids, priors)
  A <- deltas
  nd <- ncol(deltas)
  nv <- r - nd
  if (nv > 0) {
    # compute the lrlda with r - nd dimensions
    lrlda <- flashx.lrlda(X, Y, r=nv)
    A <- cbind(A, lrlda$A)
  } else if (nv < 0) {
    # retried the first r differences in the means
    A <- fm.get.cols(A, 1:r)
  }
  # otherwise nv = 0 and we exactly want the deltas
  return(list(Xr=flashx.embed(X, A), A=A))
}

flashx.rp <- function(X, r, ...) {
  X <- fm.as.matrix(X); Y <- as.vector(Y)
  n <- nrow(X); d <- ncol(X)
  A <- 1/r*fm.rnorm.matrix(d, r, 0, 1, in.mem=TRUE)
  return(list(Xr=flashx.embed(X, A), A=A))
}

flashx.lrcca <- function(X, Y, r, ...) {
  X <- fm.as.matrix(X); Y <- as.vector(Y)
  n <- length(Y); d=ncol(X)
  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }
  Yh <- lolR:::lol.utils.ohe(Y)
  S_y <- cov(Yh)
  S_yinv <- MASS::ginv(S_y)
  # covariance matrices cov(X)
  if (nrow(X) < ncol(X)) {
    Xc <- sweep(X, 2, colMeans(X), "-")
    X_cov_mul <- function(vec, extra) 1/(n-1) * t(Xc) %*% (Xc %*% vec)
    res <- fm.eigen(X_cov_mul, k=nrow(X), n=ncol(X), which="LM", sym=TRUE)
  } else {
    S_x <- cov(X)
    res <- eigen(S_x)
  }
  sigma <- res$values[res$values > 1e-6]
  U <- res$vectors[,1:length(sigma)]

  # inverse covariance matrices are ginverse in the low-rank case
  S_xy <- cov(X, fm.as.matrix(Yh))

  # decompose Sxi*Sxy*Syi*Syx
  Z <- (t(U) %*% S_xy) %*% (S_yinv %*% t(S_xy))
  Z <- t(Z)
  U <- U %*% diag(1/sigma)
  A <- tcrossprod_pca(U, Z, r)

  return(list(Xr=flashx.embed(X, A), A=A))
}

tcrossprod_pca <- function(V, U, r) {
  V <- fm.as.matrix(V)
  U <- fm.as.matrix(U)
  mul <- function(vec, extra) {
    vec <- as.vector(U %*% (t(V) %*% vec))
    vec <- V %*% (t(U) %*% vec)
  }
  res <- fm.eigen(mul, k=r, n=nrow(V), which="LM", sym=TRUE)
  res$vectors
}
