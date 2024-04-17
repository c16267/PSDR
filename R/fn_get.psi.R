#'Kernel function generator
#'@description
#'Generate a kernel function
#'@param x data matrix
#'@param y either continuous or (+1,-1) typed binary response vector
#'@param k number of basis functions for a kernel trick, floor(length(y)/3) is default
#'@return generate the kernel function \eqn{\phi(\mathbf{x})}
#' \item{\code{w}}{the first d leading eigenvectors of the matrix \eqn{\mathbf{Q}\mathbf{K}\mathbf{Q}} in Li B. et. al., (2011)}
#' \item{\code{l}}{the first d leading eigenvectors of the matrix \eqn{\mathbf{Q}\mathbf{K}\mathbf{Q}}}
#' \item{\code{scaled.x}}{scaled \eqn{\mathbf{x}}}
#' \item{\code{bw}}{bandwidth of the kernel function}
#' \item{\code{k}}{parameter for the kernel, \code{floor(length(n)/3)} is applied as a default}
#' @author Jungmin Shin, \email{jungminshin@korea.ac.kr}, Andreas Artemiou \email{artemiou@uol.ac.cy}, Seung Jun Shin, \email{sjshin@korea.ac.kr}
#'@seealso \code{\link{npsdr}}, \code{\link{phix}}
#'@import stats svmpath
#'@export get.psi

get.psi <- function(x, y, k=floor(length(y)/3)) {
  n <- nrow(x)
  x <- scale(x)
  bw <- 1/mean(as.numeric(dist(x)))^2 # bw parameter for kernel

  Kn <- svmpath::radial.kernel(x, x, bw)
  Qn <- diag(n) - matrix(1/n, n, n)

  eigen.psi <- eigen(Qn %*% Kn %*% Qn)
  Psi.new <- eigen.psi$vectors[,1:k, drop = F] # Psi
  l <- eigen.psi$values[1:k]             # eigen
  tmp.obj <- list("w"=Psi.new, "l"=l, "scaled.x"= x, "bw" = bw, "k" = k)
  #class(tmp.obj) <- "npsdr"
  return(tmp.obj)
}

