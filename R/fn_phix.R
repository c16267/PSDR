#'A function of yielding nonlinear mapping of X onto the reduced dimension
#'@description
#'A function of yielding nonlinear mapping of X onto the reduced dimension
#'@param value data matrix X
#'@param obj the object from function 'npsdr'
#'@param order the structural dimensionality
#'@return Return the value of a nonlinear Mapping X to the reduced dimension of rank 2 via the kernel function
#'@author Jungmin Shin, \email{jungminshin@korea.ac.kr}, Seung Jun Shin, \email{sjshin@korea.ac.kr}
#'@seealso \code{\link{npsdr}}
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
#'my.obj <- npsdr(x, y, H, h, lambda, delta, k=floor(length(y)/3), eps,
#'                max.iter, loss="svm")
#'x.nsvm <- phix(value=x, obj=my.obj, order=2)
#'}
#'@import stats svmpath
#'@export phix

phix <- function(value, obj, order = 2) {
  psi.function <- psi.function
  x <- value
  v <- obj$evector
  w <- obj$obj.psi$w
  l <- obj$obj.psi$l
  p <- dim(x)[2]
  kernel.function <- kernel.function(x, y=x, param.kernel = 1/p)
  tau <- mean(as.numeric(dist(x)))
  kernel.param <- 1/tau^2
  p <- ncol(x)
  if (length(value) == p) {
    temp <- psi.function(value, x, v[,1:order, drop = F], w, l, kernel.function, kernel.param)
  } else if (ncol(value) == p) {
    temp <- t(apply(value, 1, psi.function, x, v[,1:order, drop = F], w, l, kernel.function, kernel.param))
  } else if (nrow(value) == p) {
    temp <- t(apply(value, 2, psi.function, x, v[,1:order, drop = F], w, l, kernel.function, kernel.param))
  } else stop("check `str(value)`")
  temp
}


#'@noRd
#'@export
psi.function <- function(value, x, v, w, l, kernel.function, kernel.param){
  value <- matrix(value, 1, length(value))
  temp <- kernel.function(value, x, kernel.param)
  psi.value <- apply(w * c(temp - mean(temp)), 2, sum)/l
  rslt <- psi.value %*% v
  class(rslt) <- "npsdr"
  return(rslt)
}


#'@noRd
#'@export
kernel.function <- function (x, y = x, param.kernel = 1/p) {
  n <- nrow(x)
  m <- nrow(y)
  p <- ncol(x)
  normx <- drop((x^2) %*% rep(1, p))
  normy <- drop((y^2) %*% rep(1, p))
  a <- x %*% t(y)
  a <- (-2 * a + normx) + outer(rep(1, n), normy, "*")
  exp(-a * param.kernel)
}






