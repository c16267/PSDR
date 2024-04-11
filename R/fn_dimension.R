#'Determination of structural dimension
#' @description
#'Estimation of a structural dimensionality. Choose the k which maximizes a Gn value
#'@param obj The psdr object
#'@param rho hyperparameter of CVBIC
#'@return Estimated Gn scores for a optimal structural dimension will be returned
#'@author Jungmin Shin, \email{jungminshin@korea.ac.kr}, Seung Jun Shin, \email{sjshin@korea.ac.kr}
#'@seealso \code{\link{psdr}} \code{plot.Gn}
#'@examples
#'\donttest{
#'set.seed(1234)
#'n <- 200; p <- 5; H <- 10; lambda <- 0.001; eps <- 1.0e-6;
#'max.iter <- 10; h <- 1.0e-5; delta <- 2*1.0e-1;
#'init.theta <- rnorm(sd=1,n=p)
#'x <- matrix(rnorm(n*p, 0, 1), n, p)
#'err <- rnorm(n, 0, .2)
#'B <- matrix(0, p, 2)
#'B[1,1] <- 1; B[2,2] <- 1
#'x1 <- x %*% B[,1]
#'x2 <- x %*% B[,2]
#'fx <-  x1/(0.5 + (x2 + 1)^2)
#'y <- fx + err
#'obj <- psdr(x, y, init.theta, H,lambda, h, delta, eps, max.iter, loss="svm")
#'d.hat <- dimension(obj, rho=0.005); d.hat
#'plot(d.hat)
#'}
#'
#'@import stats graphics
#'@export dimension

#
# dim <- function(x, obj, rho, plot=FALSE, ...){
#     UseMethod("dim.psdr")
#}

dimension <- function(obj, rho){
  p <- nrow(obj$vectors)
  v <- obj$values
  n <- nrow(obj$x);
  temp_Mat <- diag(p); temp_Mat[lower.tri(temp_Mat)] <- 1
  temp_vec <- c(1:p)
  value <- temp_Mat %*% v - temp_vec * (rho * log(n) / n^(1/2) * v[1])
  class(value) <- "Gn"
  structure(class = "Gn", value)
  return(value)
}


#' @noRd
#' @export
plot.Gn <- function(x, ...) {
  obj <- x
  if (!inherits(obj, "Gn"))
    stop("use only with \"Gn\" objects")
  x.vec <- seq(1, length(obj), by=1)
  graphics::plot(x.vec, obj, xlab=expression(m), ylab=expression(G[n](m)), ...)
}



