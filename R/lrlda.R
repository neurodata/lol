#' Discriminant function which measures the distance between a point and a centroid j using Mahalanobis distance
#' @param x [1,r] data point that has dimensionality up to the rth subspace
#' @param centroid [1,r] centroid point j that has dimensionality up to the rth subspace
#' @param prior [1] prior probability for class j
#' @return Mahalanobis distance between x and centroid j
discriminant_fun <- function(x, centroid, prior){
  return(0.5 * sum((x - centroid)^2) - log(prior))
}

#' Linear Discriminant Analysis (LDA) Prediction
#'
#' A function for using Linear Discriminant Analysis (LDA) for prediction.
#' @param X [n1, r] data matrix with n samples in r dimensions for training.
#' @param Y [n1] the labels for each sample for training.
#' @param Xtest [n2, r] data matrix with n samples in r dimensions for testing.
#' @param Ytest [n2] the labels for each sample for testing.
#' @param ylabs [K] vector containing the unique, ordered class labels.
#' @param centroids [K, C-1] centroid matrix of the unique, ordered class labels.
#' @param priors [K] vector containing prior probability for the unique, ordered class labels.
#' @return Xr [n, K-1] projected data matrix.
#' @return Mp [K, K-1] projected centroid matrix.
#' @return Yhat [n] prediction matrix containing predictions for each example.
#' @author Richard Chen, modified by Eric Bridgeford
#' @export
fs.lda <- function(X, Y, Xtest){
  K <- length(unique(Y)) # number of classes in the data

  train_lda <- fs.project.lda(X, Y)
  ylabs <- train_lda$ylabs; priors <- train_lda$priors; cr <- train_lda$cr

  # Project the test data into the invariant subspaces
  Xr.test <- Xtest %*% train_lda$A

  # Classify the data by doing nearest centroid classification
  Yhat <- ylabs[sapply(1:(dim(Xr.test)[1]), function(i) {
    which.min(sapply(1:K, function(j) {
      discriminant_fun(Xr.test[i,], cr[j,], priors[j])
    }))
  })]

  return(Yhat)
}


#' Linear Discriminant Analysis (LDA)
#'
#' A function for implementing the LDA Algorithm.
#' @param X [n, d] the data with n samples in d dimensions.
#' @param Y [n] the labels of the samples.
#' @return A [d, r] the projection matrix
#' @return ylabs [K] vector containing the unique, ordered class labels.
#' @return centroids [K, d] centroid matrix of the unique, ordered classes.
#' @return priors [K] vector containing prior probability for the unique, ordered classes.
#' @return Xr [n, r] the data in reduced dimensionality.
#' @return cr [K, r] the centroids in reduced dimensionality.
#' @author Richard Chen, modified by Eric Bridgeford
#' @export
fs.project.lda <- function(X, Y) {
  classdat <- fs.utils.classdat(X, Y)
  priors <- classdat$priors; centroids <- classdat$centroids
  K <- classdat$K; ylabs <- classdat$ylabs
  n <- classdat$n; d <- classdat$d

  # 1. Compute the class dependent probabilities and class centroids
  priors = sapply(ylabs, function(y) sum(Y == y)/n)

  # 2. Computing with-class covariance
  W <- stats::cov(X) # within-class scatter

  # 3. Computing M* = M W^{-1/2} using the eigen-decomposition of W
  e <- eigen(W)
  V <- e$vectors # Recall that W decomposes to = V W V^T, which is V %*% diag(e$values) %*% t(V)
  W_neg_one_half <- V %*% diag(1/sqrt(e$values)) %*% t(V)
  M_star <- centroids %*% W_neg_one_half # M* = M W^{-1/2}

  # 4. Compuinge B* (which is just the covariance matrix of M*), and its eigen-decomposition, B*=V*D_BV*^T
  #    Note that the columns of V* define the coordiantes of the optimal subspaces
  B_star = stats::cov(M_star)
  V_star = eigen(B_star)$vectors[,1:K]

  # 5. Compute the full projection matrix
  A = W_neg_one_half %*% V_star

  return(list(A=A, centroids=centroids, priors=priors, ylabs=ylabs,
              Xr=X %*% A, cr=centroids %*% A))
}
