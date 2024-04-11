#'Real time sufficient dimension reduction through principal least squares SVM
#'@description
#'Real time least squares SVM SDR method for a both regression and classification problem
#'@param A a temporal matrix produced by the previous estimation with original data
#'@param r an estimated parameters from the previous estimation with original data
#'@param n a number of sample size of the original data
#'@param Xbar a mean vector from the original data
#'@param X X in new data
#'@param y y in new data, y is continuous
#'@param direction A direction of change of the new data which is either "forward" or "backward". default is "forward"
#'@param H a number of slices. default is set to 10.
#'@param lambda hyperparameter for the loss function. default is set to 0.1.
#'@return An updated psdr result will be returned in a real time fashion
#'@author Jungmin Shin, \email{jungminshin@korea.ac.kr}, Seung Jun Shin, \email{sjshin@korea.ac.kr}
#'@seealso \code{\link{psdr}}
#'@examples
#'\donttest{
#'set.seed(1234)
#'n <- 300; n2 <- 500; p <- 5; H <- 10; lambda <- 0.1;
#'x <- matrix(rnorm(n*p, 0, 1), n, p)
#'x.new <- matrix(rnorm(n2*p, 0, 1), n2, p)
#'err <- rnorm(n, 0, .2); err.new <- rnorm(n2, 0, .2)
#'B <- matrix(0, p, 2)
#'B[1,1] <- 1; B[2,2] <- 1
#'x1 <- x %*% B[,1]
#'x2 <- x %*% B[,2]
#'x1.new <- x.new %*% B[,1]
#'x2.new <- x.new %*% B[,2]
#'fx <-  x1/(0.5 + (x2 + 1)^2)
#'fx.new <-  x1.new/(0.5 + (x2.new + 1)^2)
#'y <- fx + err;                     y.new <- fx.new + err.new
#'y.binary <- c(sign(fx + err));     y.binary.new <- c(sign(fx.new + err.new))
#'ls <- psdr(x, y, H=10, lambda=0.1, loss="ls")
#'wls <- psdr(x, y.binary, H=10, lambda=0.1, loss="wls")
#'rtpsdr(A=ls$A, r=ls$r, n=ls$N, Xbar=ls$Xbar, X=x.new, y=y.new, direction="forward", H, lambda)
#'rtpsdr(A=wls$A, r=wls$r, n=wls$N, Xbar=wls$Xbar, X=x.new, y=y.binary.new,
#' direction="forward", H, lambda)
#'}
#'@import stats
#'@export rtpsdr

# -----------------------------------------------------------------------------------------------------------
# Real time Principal (weighted) least squares SDR
# -----------------------------------------------------------------------------------------------------------

rtpsdr <- function(A, r, n, Xbar, X, y, direction="forward", H=NULL, lambda=NULL){
  p <- ncol(X)
  m <- nrow(X)

  if(direction == "forward"){# data adding
    m <- nrow(X)
    if(sum(unique(y)) != 0){#plssvm
      new.Xbar <- apply(X, 2, mean)
      Xbar <- (new.Xbar*m+Xbar*n)/(n+m)
      X <- cbind(t(t(X)-Xbar), -1)
      cov.new.X <- cov(X)*(m-1)/m
      qprob <- (1:(H-1))/H
      qy <- stats::quantile(y, qprob)

      A.new <- vector(mode = "list", length = length(qprob))
      r.new <- matrix(0, ncol=p+1, nrow=length(qprob))

      if(is.null(H) == T){
        warning("H is set to 10 as a default.")
        H <- 10
      }
      if(is.null(lambda) == T){
        warning("lambda is set to 0.1 as a default.")
        lambda <- 0.1
      }
      for (s in 1:length(qprob)){
        y.tilde <- rep(1, nrow(X))
        y.tilde[y < qy[s]] <- -1  #s
        B_new <- m*cov.new.X/lambda + t(X) %*%  X
        C_new <- t(X)  %*% y.tilde
        s_part <- diag(rep(1, p+1))+A[[s]]%*%B_new
        K <- diag(rep(1, p+1))-A[[s]] %*% B_new %*% solve(s_part)

        r.new[s,] <- t(K %*% (r[s,] + A[[s]] %*% C_new))
        A.new[[s]] <- K %*% A[[s]]
      }
    }else{ #pwlssvm
      Xbar <- (apply(X, 2, mean)*m+Xbar*n)/(n+m)
      X <- cbind(t(t(X)-Xbar), -1)
      weight_list <- seq(0, 1, length=H+2)[2:(H+1)]
      A.new <- vector(mode = "list", length = H)
      r.new <- matrix(0, ncol=p+1, nrow=H)

      if(is.null(H) == T){
        warning("H is set to 10 as a default.")
        H <- 10
      }

      if(is.null(lambda) == T){
        warning("lambda is set to 0.1 as a default.")
        lambda <- 0.1
      }

      cov.new.X <- cov(X)*(m-1)/m
      for (i in 1:H){
        new.W <- diag(ifelse(y==1, 1-weight_list[i], weight_list[i]))
        B_new <- m*cov.new.X/lambda + t(X) %*% new.W %*% X
        C_new <- t(X) %*% new.W %*% y
        s_part <- diag(rep(1, p+1))+A[[i]]%*%B_new
        s <- diag(rep(1, p+1))-A[[i]] %*% B_new %*% solve(s_part)

        r.new[i,] <- t(s %*% (r[i,] + A[[i]] %*% C_new))
        A.new[[i]] <- s %*% A[[i]]
      }
    }
    Working_mat <- t(r.new[,1:p])%*%r.new[,1:p]
    eigen.Mn <- eigen(Working_mat)
    newlist <- list("x" = X[ ,-ncol(X)], "y" = y, N=n+m, Xbar=Xbar, r=r.new, A=A.new,
                    "Mn"=Working_mat, "values" = eigen.Mn$values, "vectors" = eigen.Mn$vectors)
    class(newlist) <- "psdr"
  }
  if(direction == "backward"){ #data remove
    l <- nrow(X)
    if(sum(unique(y)) != 0){#pwlssvm
      new.Xbar <- apply(X, 2, mean)
      Xbar <- (new.Xbar*l+Xbar*n)/(n+l)
      X <- cbind(t(t(X)-Xbar), -1)
      cov.new.X <- cov(X)*(l-1)/l
      qprob <- (1:(H-1))/H
      qy <- stats::quantile(y, qprob)

      A.new <- vector(mode = "list", length = length(qprob))
      r.new <- matrix(0, ncol=p+1, nrow=length(qprob))

      if(is.null(H) == T){
        warning("H is set to 10 as a default.")
        H <- 10
      }
      if(is.null(lambda) == T){
        warning("lambda is set to 0.1 as a default.")
        lambda <- 0.1
      }
      for (s in 1:length(qprob)){
        y.tilde <- rep(1, nrow(X))
        y.tilde[y < qy[s]] <- -1  #s
        B_new <- l*cov.new.X/lambda + t(X) %*%  X
        C_new <- t(X)  %*% y.tilde
        s_part <- diag(rep(1, p+1))-A[[s]]%*%B_new
        K <- diag(rep(1, p+1))+A[[s]] %*% B_new %*% solve(s_part)

        r.new[s,] <- t(K %*% (r[s,] - A[[s]] %*% C_new))
        A.new[[s]] <- K %*% A[[s]]
      }
    }else{ #pwlssvm
      Xbar <- (apply(X, 2, mean)*l+Xbar*n)/(n+l)
      X <- cbind(t(t(X)-Xbar), -1)
      weight_list <- seq(0, 1, length=H+2)[2:(H+1)]
      A.new <- vector(mode = "list", length = H)
      r.new <- matrix(0, ncol=p+1, nrow=H)

      if(is.null(H) == T){
        warning("H is set to 10 as a default.")
        H <- 10
      }

      if(is.null(lambda) == T){
        warning("lambda is set to 0.1 as a default.")
        lambda <- 0.1
      }

      cov.new.X <- cov(X)*(l-1)/l
      for (i in 1:H){
        new.W <- diag(ifelse(y==1, 1-weight_list[i], weight_list[i]))
        B_new <- l*cov.new.X/lambda + t(X) %*% new.W %*% X
        C_new <- t(X) %*% new.W %*% y
        s_part <- diag(rep(1, p+1))-A[[i]]%*%B_new
        s <- diag(rep(1, p+1))+A[[i]] %*% B_new %*% solve(s_part)

        r.new[i,] <- t(s %*% (r[i,] - A[[i]] %*% C_new))
        A.new[[i]] <- s %*% A[[i]]
      }
    }

    Working_mat <- t(r.new[,1:p])%*%r.new[,1:p]
    eigen.Mn <- eigen(Working_mat)
    newlist <- list("x" = X[ ,-ncol(X)], "y" = y, N=n+l, Xbar=Xbar, r=r.new, A=A.new,
                    "Mn"=Working_mat, "values" = eigen.Mn$values, "vectors" = eigen.Mn$vectors)
    class(newlist) <- "rtpsdr"
    structure(class = "rtpsdr", newlist)
  }
  return(newlist)
}


#' @noRd
#' @export
print.rtpsdr <- function(x, ...) {
  obj <- x
  d <- list(x = obj$x, y = obj$y, Mn= obj$Mn, evalues = obj$values, evectors = obj$vectors, N=obj$n, Xbar=  apply(obj$x, 2, mean), r=obj$r.H, A=obj$A)
  writeLines("psdr result:")
  print(d, ...)
  invisible(d)
}



#' @noRd
#' @export
plot.rtpsdr <- function(x, dim=2, ...) {
  obj <- x
  if (!inherits(obj, "psdr"))
    stop("use only with \"psdr\" objects")
  temp <- obj$vectors
  obj_psdr <- obj$x %*% temp
  par(mfrow=c(ceiling(sqrt(dim)), ceiling(sqrt(dim))))
  par(mfrow=c(1,dim))
  for(d in 1:dim){
    plot(obj_psdr[,d], obj$y, type = "p", xlab = bquote(paste(hat(B)[.(d)]^T*X)), ylab  = expression(Y) , cex=.7,...)
    graphics::lines(lowess(obj_psdr[,d], obj$y), col="red", lwd=1)
  }
  par(mfrow=c(1,1))
}
