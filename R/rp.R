#' Random Projections (RP)
#'
#' A function for implementing gaussian random projections (rp).
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param r the rank of the projection. Note that \code{r >= K}, and \code{r < d}.
#' @param scale whether to scale the random projection by the sqrt(1/d). Defaults to \code{TRUE}.
#' @param ... trailing args.
#' @return A list containing the following:
#' \item{\code{A}}{\code{[d, r]} the projection matrix from \code{d} to \code{r} dimensions.}
#' \item{\code{Xr}}{\code{[n, r]} the \code{n} data points in reduced dimensionality \code{r}.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("rp", package = "lolR")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' model <- lol.project.rp(X=X, r=5)  # use lol to project into 5 dimensions
#' @export
lol.project.rp <- function(X, r, scale=TRUE, ...){
  d <- dim(X)[2]
  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }
  # Gaussian Random Projection from d to r dimensions
  A <- array(rnorm(d*r), dim=c(d, r))

  # scale if desired
  if (scale) {
    A <- sqrt(1/r)*A
  }

  return(list(A=A, Xr=lol.embed(X, A)))
}
