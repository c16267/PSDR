#' scatter plot with the first predictor versus the response
#'@description
#' produces a scatter plot of the transformed predictors versus response variable with lowess curve
#'@param obj The gdsdr list type output
#'@param dim decide number of dimensions that will be plotted, 2 is a default
#'@param ... Generic compatibility
#'@return scatter plot with transformed mth predictor versus response
#'@author Jungmin Shin, Seung Jun Shin
#'@seealso \code{\link{psdr}}, \code{\link{dimension}}
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
#'obj <- psdr(x, y, init.theta, H,lambda, h, delta, eps, max.iter, loss="svm")
#'plot(obj, dim=2) }

#'@export
#'@import stats
#'@importFrom graphics lines plot legend


# plot.psdr <- function(x, y, ...){
#    UseMethod("plot")
# }



plot.psdr <- function(obj, dim=2, ...) {
  if (!inherits(obj, "psdr"))
    stop("use only with \"psdr\" objects")
  temp <- obj$vectors
  obj_psdr <- obj$x %*% temp
  par(mfrow=c(ceiling(sqrt(dim)), ceiling(sqrt(dim))))
  par(mfrow=c(1,dim))
  for(d in 1:dim){
    plot.default(obj_psdr[,d], obj$y, type = "p", xlab = bquote(paste(hat(B)[.(d)]^T*X)), ylab  = expression(Y) , cex=.7,...)
    graphics::lines(lowess(obj_psdr[,d], obj$y), col="red", lwd=1)
  }
  par(mfrow=c(1,1))
}


