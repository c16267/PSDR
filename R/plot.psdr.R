#' scatter plot with the first predictor versus the response
#' @description
#' produces a scatter plot of the transformed perdictors versus response variable with loess curve
#'
#'@param x X data matrix
#'@param y response variable
#'@param obj The gdsdr list type output
#'@param dim decide number of dimensions that will be plotted, 4 is a default
#'@param ... Generic compatibility
#'@return scatter plot with transformed mth predictor versus response
#'@author Jungmin Shin, Seungjun Shin
#'@seealso \code{\link{GD.psdr.solver}}, \code{\link{Gn}}
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
#'y.binary <- sign(fx + err)
#'obj <- GD.psdr.solver(x, y, init.theta, H,lambda, h, delta, eps, max.iter, loss="svm")
#'plot.psdr(x, y, obj, dim=4) }
#'@export
#'@export plot.psdr
#'@importFrom graphics lines plot legend





plot.psdr <- function(x, y, obj, dim=4, ...) {
  temp <- obj$vectors
  obj.psdr <- x %*% temp
  par(mfrow=c(ceiling(sqrt(dim)), ceiling(sqrt(dim))))
  for(d in 1:dim){
    plot(obj.psdr[,d], y, type = "p", xlab = paste("b.hat", d, "x"), ylab  = expression(Y) , cex=.5,...)
    graphics::lines(lowess(obj.psdr[,d], y), col="red", lwd=1)
  }
  par(mfrow=c(1,1))
}



