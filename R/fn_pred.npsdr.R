#'Determination of structured dimension
#' @description
#'Estimation of a structural dimensionality. Choose the k which maximizes a Gn value
#'
#'@param obj The npsdr object
#'@param new.x new data X
#'@param d true dimension for data x. d=2 is default
#'@return Estimated response variable y return
#'@author Jungmin Shin, \email{jungminshin@korea.ac.kr}, Seungjun Shin, \email{sjshin@korea.ac.kr}
#'@seealso \code{\link{GD.npsdr.solver}}
#'@examples
#'\donttest{
#'set.seed(1)
#'n <- 200;
#'p <- 5;
#'H <- 20;
#'lambda <- 0.1
#'eps <- 1.0e-5
#'max.iter <- 10
#'h <- 1.0e-6; delta <- 5*1.0e-1
#'x <- matrix(rnorm(n*p, 0, 2), n, p)
#'err <- rnorm(n, 0, .2)
#'B <- matrix(0, p, 2)
#'B[1,1] <- 1; B[2,2] <- 1
#'x1 <- x %*% B[,1]
#'x2 <- x %*% B[,2]
#'fx <-  x1/(0.5 + (x2 + 1)^2)
#'y <- c(fx + err) # response
#'y.binary <- sign(y)
#'new.x <- matrix(rnorm(n*p, 0, 2), n, p)
#'obj <- GD.npsdr.solver(x, y, H, h, lambda, delta, k=floor(length(y)/3), eps, max.iter, loss="svm")
#'pred.npsdr(obj, new.x, d = 2) }
#'@import stats MASS graphics kernlab svmpath
#'@export pred.npsdr

pred.npsdr <- function(obj, new.x, d = 2){
  x.obj <- obj$obj.psi$scaled.x
  m <- attr(x.obj, "scaled:center")
  s <- attr(x.obj, "scaled:scale")
  bw <- obj$obj.psi$bw

  w <- obj$obj.psi$w
  l <- obj$obj.psi$l
  k <- length(l)

  v <- obj$evector

  p <- length(m)
  n <- nrow(x)
  x <- matrix(x.obj, n, p)

  n.new <- length(new.x)/p

  new.x <- matrix(new.x, n.new, p)
  new.x <- t((t(new.x) - m)/s)

  Kern <- svmpath::radial.kernel(new.x, x, bw) # n * n.new
  Kern.c <- t(Kern - apply(Kern, 1, mean))
  f <- matrix(0, n.new,  k)
  for (j in 1:nrow(new.x)) {
    f[j,] <- crossprod(rep(1, n), w * Kern.c[,j])/l
  }
  pred <- f %*% v[,1:d, drop = F]
  class(pred) <- "npsdr"
  return(pred)
}



