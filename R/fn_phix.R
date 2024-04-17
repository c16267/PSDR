#'A function of yielding nonlinear mapping of \eqn{\mathbf{x}} onto the reduced dimension
#'@description
#'A function of yielding nonlinear mapping of \eqn{\mathbf{x}} onto the reduced dimension
#'@param value data matrix \eqn{\mathbf{x}}
#'@param object the object from function 'npsdr'
#'@param d the structural dimensionality. The default is \eqn{d=2}
#'@return An object with S3 class "npsdr". Details are listed below.Return the value of
#'a nonlinear Mapping \eqn{\mathbf{x}} to the reduced dimension of rank 2 via the kernel function \eqn{\phi(\mathbf{x})}
#' \item{\code{w}}{the first d leading eigenvectors of the matrix \eqn{\mathbf{Q}\mathbf{K}\mathbf{Q}} in Li B. et. al., (2011)}
#' \item{\code{l}}{the first d leading eigenvectors of the matrix \eqn{\mathbf{Q}\mathbf{K}\mathbf{Q}}}
#' \item{\code{scaled.x}}{scaled \eqn{\mathbf{x}}}
#' \item{\code{bw}}{bandwidth of the kernel function}
#' \item{\code{k}}{parameter for the kernel, \code{floor(length(n)/3)} is applied as a default}
#'@references Li, B., Artemiou, A. and Li, L. (2011)
#' \emph{Principal support vector machines for linear and
#' nonlinear sufficient dimension reduction, Annals of Statistics 39(6): 3182–3210}.
#' @author Jungmin Shin, \email{jungminshin@korea.ac.kr}, Andreas Artemiou \email{artemiou@uol.ac.cy}, Seung Jun Shin, \email{sjshin@korea.ac.kr}
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
#'obj <- npsdr(x, y, H, h, lambda, delta, k=floor(length(y)/3), eps,
#'                max.iter, loss="svm")
#'phix(value=x, object=obj, d=2)
#'}
#'@import stats svmpath
#'@export phix

phix <- function(value, object, d = 2) {
  psi.function <- psi.function
  x <- value
  v <- object$evector
  w <- object$obj.psi$w
  l <- object$obj.psi$l
  p <- dim(x)[2]
  kernel.function <- kernel.function(x, y=x, param.kernel = 1/p)
  tau <- mean(as.numeric(dist(x)))
  kernel.param <- 1/tau^2
  p <- ncol(x)
  if (length(value) == p) {
    temp <- psi.function(value, x, v[,1:d, drop = F], w, l, kernel.function, kernel.param)
  } else if (ncol(value) == p) {
    temp <- t(apply(value, 1, psi.function, x, v[,1:d, drop = F], w, l, kernel.function, kernel.param))
  } else if (nrow(value) == p) {
    temp <- t(apply(value, 2, psi.function, x, v[,1:d, drop = F], w, l, kernel.function, kernel.param))
  } else stop("check `str(value)`")
  temp
}


#'@noRd
psi.function <- function(value, x, v, w, l, kernel.function, kernel.param){
  value <- matrix(value, 1, length(value))
  temp <- kernel.function(value, x, kernel.param)
  psi.value <- apply(w * c(temp - mean(temp)), 2, sum)/l
  rslt <- psi.value %*% v
  class(rslt) <- "npsdr"
  return(rslt)
}


#'@noRd
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






