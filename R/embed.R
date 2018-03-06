#' Embedding
#'
#' A function that embeds points in high dimensions to a lower dimensionality.
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param A \code{[d, r]} the embedding matrix from \code{d} to \code{r} dimensions.
#' @param ... optional args.
#' @return an array \code{[n, r]} the original \code{n} points embedded into \code{r} dimensions.
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.lol(X=X, Y=Y, r=5)  # use lol to project into 5 dimensions
#' Xr <- lol.embed(X, model$A)
#' @author Eric Bridgeford
#' @export
lol.embed <- function(X, A, ...) {
  return(X %*% A)
}
