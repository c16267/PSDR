#' Unified Principal sufficient dimension reduction methods
#'
#' A unified and user-friendly \proglang{R} package for applying the principal sufficient dimension reduction methods for both linear and nonlinear ,and regression and classification context.
#' The package has an extendable power by varying loss functions for the SVM, even for an user-defined arbitrary function,
#' unless those are convex and differentiable everywhere over the support.
#' Also, it provides a realtime sufficient dimension reduction update procedure using the principal least squares SVM.
#'
#' Details on \code{loss} option:
#'
#' The argument \code{loss} determines a specific loss function for SVM and the corresponding SDR method. For example, \code{loss="ls"} means that the user can do SDR with
#' least square SVM. The package provides several pre-embeded functions. 1. for regression problem: \code{loss="svm"} is the hinge loss, \code{loss="logstic"} is the logistic loss, \code{loss="l2svm"} is the squared hinge loss,
#' \code{loss="LUM"} is for the large margin unified loss, \code{loss="asymls"} is for asymmetric least square loss
#' 2. Also the corresponding weighted loss functions are included, such as, \code{loss="wsvm"} , \code{loss="wlogistic"}, \code{loss="l2wsvm"},
#' \code{loss="wLUM"} and \code{loss="wls"}, which mean weighted hinge loss, weighted logistic loss, weighted squared hinge loss, weighted LUM loss and weighted least square loss, respectively.
#'
#' Not only function \code{psdr} includes popular loss functions, but also, it is designed for working with user defined arbitrary convex loss function that is claimed through the argument \code{loss}.
#' Two examples of the usage of user-defined losses are presented below (\code{m} represents a margin):
#'
#' \code{myLogistic <- function(m,...){rslt <- log(1+exp(-m)) return(rslt)}},
#'
#' \code{myLS <- function(m,...){rslt <- (1-m)^2 return(rslt)}}.
#'
#' Users can define their own loss function in advance, and apply those functions to \code{psdr} like \code{loss="myLogistic"}
#' or \code{loss="myLS"}.
#'
#' @param x input matrix, of dimension \code{nobs} x \code{nvars}; each row is an observation vector. Requirement: \code{nvars}>1; in other words, \code{x} should have 2 or more columns.
#' @param y response variable, either can be continuous variable or (+1,-1) coded binary response vector
#' @param init initial coefficient vector of which dimension matches the \code{nvar} of \eqn{\mathbf{X}}. If it is not specified, random vector from standard normal distribution is applied by default
#' @param H the number of slices and probabilities equally spaced in \eqn{(0,1)}. default value is 10
#' @param lambda the cost parameter for the svm loss function. The default value is 0.1
#' @param delta learning rate for the gradient descent algorithm. The default value is 0.1
#' @param h very small interval for calculating numerical derivatives for a given arbitrary loss function. The default is 1.0e-5
#' @param eps the threshold for stopping iteration with respect to the magnitude of the change of the derivative. The default value is 1.0e-5
#' @param max.iter maximum iteration number for the optimization process. default value is 30
#' @param loss pre-specified loss functions are "logistic", "svm","l2svm","wsvm", etc. and user-defined loss function object also can be used formed by inside double quotation mark
#' @param a the first hyperparameter for the LUM loss function
#' @param c the second hyperparameter for the LUM loss function
#' @param stochastic If \code{TRUE} then the stochastic gradient descent algorithm will be implemented to optimize the loss function. The default is FALSE
#' @param plot If \code{TRUE} then it produces scatter plots of \eqn{Y} versus \eqn{\hat{B^{\top}}_{j}\mathbf{X}}. \eqn{j} can be specified by the user with \eqn{j=2} as a default. The default is FALSE
#'
#' @return An object with S3 class "psdr". Details are listed below.
#' \item{\code{Mn}}{The estimated working matrix, which is obtained by the cumulative
#' outer product of the estimated parameters over H}
#' \item{\code{evalues}}{Eigenvalues of the Mn}
#' \item{\code{evectors}}{Eigenvectors of the Mn, the first leading d eigenvectors consists
#' the basis of the central subspace}
#' @author Jungmin Shin, \email{jungminshin@korea.ac.kr}, Andreas Artemiou \email{artemiou@uol.ac.cy}, Seung Jun Shin, \email{sjshin@korea.ac.kr}
#' @references Artemiou, A. and Dong, Y. (2016)
#' \emph{Sufficient dimension reduction via principal lq support vector machine,
#'  Electronic Journal of Statistics 10: 783–805}.\cr
#'  Artemiou, A., Dong, Y. and Shin, S. J. (2021)
#' \emph{Real-time sufficient dimension reduction through principal least
#'  squares support vector machines, Pattern Recognition 112: 107768}.\cr
#'  Kim, B. and Shin, S. J. (2019)
#' \emph{Principal weighted logistic regression for sufficient dimension
#' reduction in binary classification, Journal of the Korean Statistical Society 48(2): 194–206}.\cr
#'  Li, B., Artemiou, A. and Li, L. (2011)
#' \emph{Principal support vector machines for linear and
#' nonlinear sufficient dimension reduction, Annals of Statistics 39(6): 3182–3210}.\cr
#' Soale, A.-N. and Dong, Y. (2022)
#' \emph{On sufficient dimension reduction via principal asymmetric
#'  least squares, Journal of Nonparametric Statistics 34(1): 77–94}.\cr
#'  Wang, C., Shin, S. J. and Wu, Y. (2018)
#' \emph{Principal quantile regression for sufficient dimension
#'  reduction with heteroscedasticity, Electronic Journal of Statistics 12(2): 2114–2140}.\cr
#'  Shin, S. J., Wu, Y., Zhang, H. H. and Liu, Y. (2017)
#' \emph{Principal weighted support vector machines for sufficient dimension reduction in
#'  binary classification, Biometrika 104(1): 67–81}. \cr
#'  Li, L. (2007)
#' \emph{Sparse sufficient dimension reduction, Biometrika 94(3): 603–613}.
#' @seealso \code{\link{plot}}, \code{\link{print}}, \code{\link{crBIC}}
#' @examples
#' set.seed(1)
#' n <- 200;
#' p <- 5;
#' H <- 10;
#' lambda <- 0.1
#' eps <- 1.0e-5
#' max.iter <- 30
#' init.theta <- rnorm(p,0,1)
#' h <- 1.0e-5; delta <- 0.5
#' x <- matrix(rnorm(n*p, 0, 2), n, p)
#' err <- rnorm(n, 0, .2)
#' B <- matrix(0, p, 2)
#' B[1,1] <- 1; B[2,2] <- 1
#' x1 <- x %*% B[,1]
#' x2 <- x %*% B[,2]
#' fx <-  x1/(0.5 + (x2 + 1)^2)
#' y <- c(fx + err) # response
#' my.hinge <- function(m,...){
#'   rslt <- (1-m)*(as.numeric((1-m) > 0))
#'   return(rslt)
#' }
#' obj <- psdr(x, y, init.theta, H, lambda, delta, h, eps, max.iter, loss="svm")
#' psdr(x, y, init.theta, H,lambda, delta, h, eps, max.iter, loss="my.hinge")
#' print(obj)
#' plot(obj)
#'
#' ##real data: Boston housing data
#' data("BostonHousing")
#' attach(BostonHousing)
#' BostonHousing <- BostonHousing[BostonHousing$crim < 3.2 , -c(4,9)]
#' X <- BostonHousing[,-12]
#' Y <- BostonHousing[,"medv"]
#' p <- ncol(X); H <- 20; lambda <- 0.1; eps <- 1.0e-5;
#' max.iter <- 100; h <- 1.0e-5; delta <- 2*1.0e-1;
#' set.seed(1); init.theta <- rnorm(sd=1,n=p)
#' rslt <- psdr(X, Y, init.theta, H,lambda, h, delta, eps, max.iter, loss="svm")
#' value.lsvm <- rslt$values
#' lsvm <- round(rslt$vectors,3)
#' X <- as.matrix(X)
#' x.lsvm <- X %*% lsvm
#' plot(x.lsvm[,1],Y , type = "p", xlab = expression(hat(b)[1]^T*X), ylab="medv", cex=1)
#' lines(lowess( x.lsvm[,1], Y), col="red", lwd=2)
#' plot(x.lsvm[,2], Y, type = "p", xlab = expression(hat(b)[2]^T*X), ylab="medv", cex=1);
#' lines(lowess(x.lsvm[,2], Y), col="blue", lwd=2)
#'
#' ##real data: Wisconsin diagnostic breast cancer data
#' data(wisc)
#' x.wisc <- matrix(unlist(wisc[,-c(1,2)]), ncol = 30)
#' y.wisc <- 2*as.numeric(as.factor(unlist(wisc[,2]))) - 3
#' init.theta <- rnorm(dim(x.wisc)[2],0,1)
#' wisc.obj <- psdr(x.wisc, y.wisc, init.theta, H=20,lambda=0.1, h=1.0e-6,
#'                  delta=0.5,eps=10^-4, max.iter=30, loss="wlogistic")
#' value.lsvm <- wisc.obj$values
#' lsvm <- round(wisc.obj$vectors,3)
#' x.lsvm <- x.wisc %*% lsvm
#' par(mar=c(5,5,5,5), oma=c(1,1,1,1))
#' plot(x.lsvm[,1], x.lsvm[,2], type = "n", xlab = expression(hat(b)[1]^T*X), ylab = expression(hat(b)[2]^T*X))
#' points(x.lsvm[y.wisc == 1,1], x.lsvm[y.wisc == 1,2], col = 4, pch = "+")
#' points(x.lsvm[y.wisc != 1,1], x.lsvm[y.wisc != 1,2], col = 2)
#'@import stats
#'@export psdr



psdr <- function(x, y, init=NULL, H=NULL, lambda=NULL, delta=NULL, h=1.0e-5, eps=1.0e-5, max.iter=NULL,
                 loss=NULL, a=NULL, c=NULL, stochastic=FALSE, plot=FALSE) {
  if(sum(as.character(loss) == c("ls", "wls")) == 0){
    if(!is.matrix(x) & !is.data.frame(x))
      stop("x must be a matrix or dataframe.")
    if(ncol(as.matrix(y)) != 1)
      stop("y must be a univariate.")
    if(is.null(delta) == T ){
      writeLines("Delta is a learning rate which should be specified as a positve value. It is set to 0.1")
      delta <- 0.1
    }

    if(is.null(loss)==T){
      stop("Loss function should be called.")
      #loss <- 'svm'
    }

    if(as.character(loss) == "LUM" ){
      if(is.null(a) == T)
        warning("A LUM loss fucntion needs both hyperparameter a and c.")
      #a <- 1; c<- 10^8
    }

    if(is.null(H) == T){
      writeLines("H is set to 10 as a default.")
      H <- 10
    }
    if(is.null(h) == T){
      writeLines("h is not specifed. It is set to 1.0e-5 as a default.")
      h <- 1.0e-5
    }
    if(is.null(max.iter) == T){
      writeLines("max.iter is set to 30 as a default.")
      max.iter <- 30
    }
    if(is.null(init) == T){
      init <- rnorm(sd=1,n=p)
      writeLines("initial parameter is generated from N(0,1).")
    }
    if(length(init) != ncol(x))
      stop("a dimension of the initial theta must match that of input x.")
  }


  if(is.null(lambda) == T){
    writeLines("lambda is set to 1 as a default.")
    lambda <- 1
  }
  if (length(y) != nrow(x)){     #check lengths
    stop("The response and predictors have different number of observations.")}

  n <- nrow(x)
  p <- ncol(x)

  bar.x <- apply(x, 2, mean)
  x.star <- cbind(t(t(x)-bar.x), -1)
  cov.x.star <- cov(x.star)
  cov.x <- cov(x)
  step <- 1/H
  pi.grid <- seq(step, 1-step, by = step)

  # generate y.tilde
  qprob <- (1:(H-1))/H
  qy <- stats::quantile(y, qprob)

  # standardization
  temp <- eigen(cov.x)
  D <- diag(sqrt(temp$values))
  V <- temp$vectors
  sd.x <-  V %*% D
  inv.sd.x <- diag(1/(sqrt(temp$values))) %*% t(V)

  centered.x <- t(x) - bar.x  # p x n
  z.new <- t((inv.sd.x %*% centered.x))
  eigen.mat <- diag(1,p,p)
  theta.new <- rep(0,p)
  type.list <- c("svm","logistic","l2svm", "wsvm", "LUM", "quantile","asymls", "wlogistic","l2wsvm","wLUM")
  type.list2 <- c("ls","wls")

  if(sum(as.character(loss) == type.list) != 0){
    w.init <- matrix(init, nrow=p, ncol=length(qprob))
    w.final <- matrix(0, nrow=p, ncol=length(qprob))

    if(as.character(loss) == "svm"){
      if(sum(unique(y)) == 0){
        stop("response variable should be continuous!")}

      for (s in 1:length(qprob)) {
        y.tilde.new <- rep(1, nrow(x))
        y.tilde.new[y < qy[s]] <- -1  #s
        pos.rate <- sum(y.tilde.new==1)/nrow(z.new)
        neg.rate <- sum(y.tilde.new==-1)/nrow(z.new)
        for(iter in 1:max.iter){
          if(stochastic == TRUE){  ###stratified sampling
            set.seed(iter)
            pos.ind <- sample(which(y.tilde.new==1), ceiling((nrow(z.new)/(log(nrow(z.new))))*pos.rate))
            neg.ind <- sample(which(y.tilde.new==-1), floor((nrow(z.new)/(log(nrow(z.new))))*neg.rate))
            ind <- c(pos.ind,neg.ind)
            z <- z.new[ind,]
            y.tilde <- y.tilde.new[ind]
          }else{
            z <- z.new
            y.tilde <- y.tilde.new
          }
          n <- nrow(z)
          w <- w.init
            for (k in 1:p){
              margin.v <- z %*% w[,s] * y.tilde #s
              deriv <- -z[,k]*y.tilde*as.numeric(I((1-margin.v)>0)) #k
              derivative.j <- lambda*mean(deriv) + 2*w[k,s]  ##k,s
              theta.new[k] <- w[k,s] -  delta*derivative.j  ##k,k,s
            }
            w[,s] <- theta.new
            w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
            if(max(abs(deriv)) < eps)
              break
          }
          w.final[,s] <- w[,s]
        }
      }
    if(as.character(loss) == "l2svm"){
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
      for (s in 1:length(qprob)) {
        y.tilde.new <- rep(1, nrow(x))
        y.tilde.new[y < qy[s]] <- -1  #s
        pos.rate <- sum(y.tilde.new==1)/nrow(z.new)
        neg.rate <- sum(y.tilde.new==-1)/nrow(z.new)
        for(iter in 1:max.iter){
          if(stochastic == TRUE){  ###stratified sampling
            set.seed(iter)
            pos.ind <- sample(which(y.tilde.new==1), ceiling((nrow(z.new)/(log(nrow(z.new))))*pos.rate))
            neg.ind <- sample(which(y.tilde.new==-1), floor((nrow(z.new)/(log(nrow(z.new))))*neg.rate))
            ind <- c(pos.ind,neg.ind)
            z <- z.new[ind,]
            y.tilde <- y.tilde.new[ind]
          }else{
            z <- z.new
            y.tilde <- y.tilde.new
          }
          n <- nrow(z)
          w <- w.init
          for (k in 1:p){
            margin.v <- z %*% w[,s] * y.tilde #s
            deriv <- -z[,k]*y.tilde*as.numeric(I((1-margin.v)>0))*2*(1-margin.v) #k
            derivative.j <- lambda*mean(deriv) + 2*w[k,s]  ##k,s
            theta.new[k] <- w[k,s] -  delta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new
          w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
          if(max(abs(deriv)) < eps)
            break
        }
        w.final[,s] <- w[,s]
      }
    }
    if(as.character(loss) == "l2wsvm"){
      if(sum(unique(y)) != 0)
        stop("response variable should be a binary type!")
      y.new <- y
      pos.rate <- sum(y.new==1)/nrow(z.new)
      neg.rate <- sum(y.new==-1)/nrow(z.new)
      for (s in 1:length(qprob)) {
        for(iter in 1:max.iter){
          if(stochastic == TRUE){  ###stratified sampling
            set.seed(iter)
            pos.ind <- sample(which(y.new==1), ceiling((nrow(z.new)/(log(nrow(z.new))))*pos.rate))
            neg.ind <- sample(which(y.new==-1), floor((nrow(z.new)/(log(nrow(z.new))))*neg.rate))
            ind <- c(pos.ind,neg.ind)
            z <- z.new[ind,]
            y.bi <- y.new[ind]
          }else{
            z <- z.new
            y.bi <- y.new
          }
          n <- nrow(z)
          w <- w.init
          for (k in 1:p){
            margin.v <- z %*% w[,s] * y.bi #s
            weight <- (1-pi.grid[s])*as.numeric(y.bi==1)+(pi.grid[s])*(as.numeric(y.bi==-1)) #s,s
            deriv <- -z[,k]*y.bi*as.numeric(I((1-margin.v)>0))*2*(1-margin.v)*weight #k
            derivative.j <- lambda*mean(deriv) + 2*w[k,s]  ##k,s
            theta.new[k] <- w[k,s] -  delta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new  #s
          w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
          if(max(abs(deriv)) < eps)
            break
        }
        w.final[,s] <- w[,s]
      }
    }

    if(as.character(loss) == "logistic"){
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
      for (s in 1:length(qprob)) {
        y.tilde.new <- rep(1, nrow(x))
        y.tilde.new[y < qy[s]] <- -1  #s
        pos.rate <- sum(y.tilde.new==1)/nrow(z.new)
        neg.rate <- sum(y.tilde.new==-1)/nrow(z.new)
        for(iter in 1:max.iter){
          if(stochastic == TRUE){  ###stratified sampling
            set.seed(iter)
            pos.ind <- sample(which(y.tilde.new==1), ceiling((nrow(z.new)/(log(nrow(z.new))))*pos.rate))
            neg.ind <- sample(which(y.tilde.new==-1), floor((nrow(z.new)/(log(nrow(z.new))))*neg.rate))
            ind <- c(pos.ind,neg.ind)
            z <- z.new[ind,]
            y.tilde <- y.tilde.new[ind]
          }else{
            z <- z.new
            y.tilde <- y.tilde.new
          }
          n <- nrow(z)
          w <- w.init
          for (k in 1:p){
            margin.v <- z %*% w[,s] * y.tilde #s
            deriv <- -z[,k]*y.tilde*(1/(1+exp(margin.v)))#k
            derivative.j <- lambda*mean(deriv) + 2*w[k,s]  ##k,s
            theta.new[k] <- w[k,s] -  delta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new
          w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
          if(max(abs(deriv)) < eps)
            break
        }
        w.final[,s] <- w[,s]
      }
    }

    if(as.character(loss) == "wlogistic"){
      for (s in 1:length(qprob)){
        if(sum(unique(y)) != 0)
          stop("response variable should be a binary type!")
        y.new <- y
        pos.rate <- sum(y.new==1)/nrow(z.new)
        neg.rate <- sum(y.new==-1)/nrow(z.new)
        for(iter in 1:max.iter){
          if(stochastic == TRUE){  ###stratified sampling
            set.seed(iter)
            pos.ind <- sample(which(y.new==1), ceiling((nrow(z.new)/(log(nrow(z.new))))*pos.rate))
            neg.ind <- sample(which(y.new==-1), floor((nrow(z.new)/(log(nrow(z.new))))*neg.rate))
            ind <- c(pos.ind,neg.ind)
            z <- z.new[ind,]
            y.bi <- y.new[ind]
          }else{
            z <- z.new
            y.bi <- y.new
          }
          n <- nrow(z)
          w <- w.init
          for (k in 1:p){
            margin.v <- z %*% w[,s] * y.bi #s
            weight <- (1-pi.grid[s])*as.numeric(y.bi==1)+(pi.grid[s])*(as.numeric(y.bi==-1)) #s,s
            deriv <- weight*(-z[,k])*y.bi*(1/(1+exp(margin.v)))   #k
            derivative.j <- lambda*mean(deriv) + 2*w[k,s]  ##k,s
            theta.new[k] <- w[k,s] -  delta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new
          w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
          if(max(abs(deriv)) < eps)
            break
        }
        w.final[,s] <- w[,s]
      }
    }

    if(as.character(loss) == "wsvm"){
      if(sum(unique(y)) != 0)
        stop("response variable should be a binary type!")
      y.new <- y
      pos.rate <- sum(y.new==1)/nrow(z.new)
      neg.rate <- sum(y.new==-1)/nrow(z.new)
      for (s in 1:length(qprob)) {
        for(iter in 1:max.iter){
          if(stochastic == TRUE){  ###stratified sampling
            set.seed(iter)
            pos.ind <- sample(which(y.new==1), ceiling((nrow(z.new)/(log(nrow(z.new))))*pos.rate))
            neg.ind <- sample(which(y.new==-1), floor((nrow(z.new)/(log(nrow(z.new))))*neg.rate))
            ind <- c(pos.ind,neg.ind)
            z <- z.new[ind,]
            y.bi <- y.new[ind]
          }else{
            z <- z.new
            y.bi <- y.new
          }
          n <- nrow(z)
          w <- w.init
          for (k in 1:p){
            margin.v <- z %*% w[,s] * y.bi #s
            weight <- (1-pi.grid[s])*as.numeric(y.bi==1)+(pi.grid[s])*(as.numeric(y.bi==-1)) #s,s
            deriv <- -weight*z[,k]*y.bi*as.numeric(I((1-margin.v)>0)) #k
            derivative.j <- lambda*mean(deriv) + 2*w[k,s]  ##k,s
            theta.new[k] <- w[k,s] -  delta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new  #s
          w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
          if(max(abs(deriv)) < eps)
            break
        }
        w.final[,s] <- w[,s]
      }
    }

    if(as.character(loss) == "LUM"){
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
      for (s in 1:length(qprob)) {
        y.tilde.new <- rep(1, nrow(x))
        y.tilde.new[y < qy[s]] <- -1  #s
        pos.rate <- sum(y.tilde.new==1)/nrow(z.new)
        neg.rate <- sum(y.tilde.new==-1)/nrow(z.new)
        for(iter in 1:max.iter){
          if(stochastic == TRUE){  ###stratified sampling
            set.seed(iter)
            pos.ind <- sample(which(y.tilde.new==1), ceiling((nrow(z.new)/(log(nrow(z.new))))*pos.rate))
            neg.ind <- sample(which(y.tilde.new==-1), floor((nrow(z.new)/(log(nrow(z.new))))*neg.rate))
            ind <- c(pos.ind,neg.ind)
            z <- z.new[ind,]
            y.tilde <- y.tilde.new[ind]
          }else{
            z <- z.new
            y.tilde <- y.tilde.new
          }
          n <- nrow(z)
          w <- w.init

          for (k in 1:p){
            margin.v <- z %*% w[,s] * y.tilde #s
            deriv <- -z[,k]*y.tilde*as.numeric(margin.v < c/(1+c)) +
              (a/(1+c))*{(a/((1+c)*margin.v-c+a))^(a-1)}*( (1+c)*a*(-z[,k]*y.tilde)) / {(1+c)*margin.v-c+a}^2*
              as.numeric(margin.v >= c/(1+c))
            derivative.j <- lambda*mean(deriv) + 2*w[k,s]  ##k,s
            theta.new[k] <- w[k,s] -  delta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new
          w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
          if(max(abs(deriv)) < eps)
            break
        }
        w.final[,s] <- w[,s]
      }
    }
    if(as.character(loss) == "wLUM"){
      if(sum(unique(y)) != 0)
        stop("response variable should be a binary type!")
      y.new <- y
      pos.rate <- sum(y.new==1)/nrow(z.new)
      neg.rate <- sum(y.new==-1)/nrow(z.new)
      for (s in 1:length(qprob)) {
        for(iter in 1:max.iter){
          if(stochastic == TRUE){  ###stratified sampling
            set.seed(iter)
            pos.ind <- sample(which(y.new==1), ceiling((nrow(z.new)/(log(nrow(z.new))))*pos.rate))
            neg.ind <- sample(which(y.new==-1), floor((nrow(z.new)/(log(nrow(z.new))))*neg.rate))
            ind <- c(pos.ind,neg.ind)
            z <- z.new[ind,]
            y.bi <- y.new[ind]
          }else{
            z <- z.new
            y.bi <- y.new
          }
          n <- nrow(z)
          w <- w.init
          for (k in 1:p){
            margin.v <- z %*% w[,s] * y.bi #s
            weight <- (1-pi.grid[s])*as.numeric(y.bi==1)+(pi.grid[s])*(as.numeric(y.bi==-1)) #s,s
            deriv <- -z[,k]*y.bi*as.numeric(margin.v < c/(1+c)) +
              (a/(1+c))*{(a/((1+c)*margin.v-c+a))^(a-1)}*( (1+c)*a*(-z[,k]*y.bi)) / {(1+c)*margin.v-c+a}^2*
              as.numeric(margin.v >= c/(1+c))*weight
            derivative.j <- lambda*mean(deriv) + 2*w[k,s]  ##k,s
            theta.new[k] <- w[k,s] -  delta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new  #s
          w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
          if(max(abs(deriv)) < eps)
            break
        }
        w.final[,s] <- w[,s]
      }
    }
    if(as.character(loss) == "quantile"){
      #if(sum(unique(y)) == 0)
      #  stop("response variable should be continuous!")
      for (s in 1:length(pi.grid)) {
        for(iter in 1:max.iter){
          if(stochastic == TRUE){  ###stratified sampling
            set.seed(iter)
            ind <- sample(length(y), sqrt(p*length(y))) #sqrt(np)
            z <- z.new[ind,]
            y.new <- y[ind]
          }else{
            z <- z.new
            y.new <- y
          }
          n <- nrow(z)
          w <- w.init
          for (k in 1:p){
            u <- y.new - z%*%w[,s]   #s
            #derivative.j <- 2*w[k,s]+lambda*sum((-z[,k]*pi.grid[s])*as.numeric(I(u>0))+(z[,k]*(1-pi.grid[s]))*as.numeric(I(u<0)))*(1/length(y.new)) #k,s,k,s,k,s
            derivative.j <- 2*w[k,s]+lambda*(1/length(y))*sum(-z[,k]*{pi.grid[s]*as.numeric(I(u>0))+(1-pi.grid[s])*as.numeric(I(u<=0))})
            theta.new[k] <- w[k,s] -  delta*derivative.j  ##k,k,s
            #if(derivative.j < eps)
            #  break
          }
          #print(derivative.j)
          w[,s] <- theta.new  #s
          w.init <- matrix(theta.new, nrow=p, ncol = length(pi.grid))
        }
        w.final[,s] <- w[,s] #s,s
      }
    }
    if(as.character(loss) == "asymls"){
      # if(sum(unique(y)) == 0)
      #   stop("response variable should be continuous!")
      for (s in 1:length(pi.grid)) {
        for(iter in 1:max.iter){
          if(stochastic == TRUE){  ###stratified sampling
            set.seed(iter)
            ind <- sample(length(y), sqrt(p*length(y))) #sqrt(np)
            z <- z.new[ind,]
            y.new <- y[ind]
          }else{
            z <- z.new
            y.new <- y
          }
          n <- nrow(z)
          w <- w.init
          for (k in 1:p){
            u <- y.new - z%*%w[,s]
            derivative.j <- 2*w[k,s]+lambda*(1/length(y))*sum((-z[,k]*2*u)*{pi.grid[s]*as.numeric(I(u>=0))+(1-pi.grid[s])*as.numeric(I(u<=0))})
            theta.new[k] <- w[k,s] -  delta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new  #s
          w.init <- matrix(theta.new, nrow=p, ncol = length(pi.grid))
        }
        w.final[,s] <- w[,s] #s,s
      }
    }

    psi <- t(inv.sd.x) %*% w.final
    Mn <- matrix(0, p, p)
    for (h in 1:length(qprob)) Mn <- Mn + psi[,h, drop = F] %*% t(psi[,h, drop = F])
    eigen.Mn <- eigen(Mn)
    newlist <- list("x" = x, "y" = y,"Mn"=Mn, "values" = eigen.Mn$values, "vectors" = eigen.Mn$vectors)
    structure(class = "psdr", newlist)
    class(newlist) <- "psdr"
    return(newlist)
   }

  else if(sum(as.character(loss) == type.list) == 0 & sum(as.character(loss)==type.list2)==0){
    if(h > 1.0e-3)
      writeLines("An infinitesimal interval h should be smaller than 1.0e-3.")
    ft <- E(loss)
    grid.m <- seq(-2,2,length=100)
    if(plot == TRUE){
      plot(grid.m, ft(grid.m,a,c,prob=0.5), type="l", xlab="margin", ylab="loss")
    }
    writeLines("loss function must be a convex function")
    w.init <- matrix(init, nrow=p, ncol=length(qprob))
    w.final <- matrix(0, nrow=p, ncol=length(qprob))

    ## arbitrary loss continuous ##
    if(sum(unique(y)) != 0){
      for(s in 1:length(qprob)){
        y.tilde.new <- rep(1, nrow(x))
        y.tilde.new[y < qy[s]] <- -1  #s
        pos.rate <- sum(y.tilde.new==1)/nrow(z.new)
        neg.rate <- sum(y.tilde.new==-1)/nrow(z.new)
        w <- w.init
        for(iter in 1:max.iter){
          if(stochastic == TRUE){  ###stratified sampling
            set.seed(iter)
            pos.ind <- sample(which(y.tilde.new==1), ceiling((nrow(z.new)/(log(nrow(z.new))))*pos.rate))
            neg.ind <- sample(which(y.tilde.new==-1), floor((nrow(z.new)/(log(nrow(z.new))))*neg.rate))
            ind <- c(pos.ind,neg.ind)
            z <- z.new[ind,]
            y.tilde <- y.tilde.new[ind]
          }else{
            z <- z.new
            y.tilde <- y.tilde.new
          }
          n <- nrow(z)
          derivative.vec <- rep(0,p)
          for (k in 1:p){
            star <- (fn_arbitrary_loss(z, y.tilde,theta=w[,s]+h*eigen.mat[,k],lambda,loss,a,c,prob=pi.grid[s])-
                       fn_arbitrary_loss(z, y.tilde,theta=w[,s]-h*eigen.mat[,k],lambda,loss,a,c,prob=pi.grid[s]))
            derivative.vec[k] <- sign(star)*exp( log(abs(star))-log(2*h))
            theta.new[k] <- w[k,s] - delta * derivative.vec[k] ###k,k,s,k
          }
          w[,s] <- theta.new   #s
          if(max(abs(derivative.vec)) < eps)
            break
        }
        w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
        w.final[,s] <- w[,s] ##s,s
      }
    }

    ## arbitrary loss binary##
    if(sum(unique(y))==0){
      y.new <- y
      pos.rate <- sum(y.new==1)/nrow(z.new)
      neg.rate <- sum(y.new==-1)/nrow(z.new)
      for(s in 1:length(qprob)){
        w <- w.init
        for(iter in 1:max.iter){
          if(stochastic == TRUE){  ###stratified sampling
            set.seed(iter)
            pos.ind <- sample(which(y.new==1), ceiling((nrow(z.new)/(log(nrow(z.new))))*pos.rate))
            neg.ind <- sample(which(y.new==-1), floor((nrow(z.new)/(log(nrow(z.new))))*neg.rate))
            ind <- c(pos.ind,neg.ind)
            z <- z.new[ind,]
            y.bi <- y.new[ind]
          }else{
            z <- z.new
            y.bi <- y.new
          }
          derivative.vec <- rep(0,p)
          for (k in 1:p){
            star <- (fn_arbitrary_binary_loss(z, y.bi, prob=pi.grid[s],theta=w[,s]+h*eigen.mat[,k],lambda,loss,a,c)-
                       fn_arbitrary_binary_loss(z, y.bi,prob=pi.grid[s],theta=w[,s]-h*eigen.mat[,k],lambda,loss,a,c))
            derivative.vec[k] <- sign(star)*exp( log(abs(star))-log(2*h))
            theta.new[k] <- w[k,s] - delta * derivative.vec[k]
          }
          w[,s] <- theta.new
          if(max(abs(derivative.vec)) < eps)
            break
        }
        w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
        w.final[,s] <- w[,s]
      }
    }
    psi <- t(inv.sd.x) %*% w.final
    Mn <- matrix(0, p, p)
    for (h in 1:length(qprob)) Mn <- Mn + psi[,h, drop = F] %*% t(psi[,h, drop = F])
    eigen.Mn <- eigen(Mn)
    newlist <- list("x" = x, "y" = y, "Mn"=Mn, "values" = eigen.Mn$values, "vectors" = eigen.Mn$vectors)

    class(newlist) <- "psdr"
    structure(class = "psdr", newlist)
    return(newlist)
  }
  else if(sum(as.character(loss) == type.list2) != 0){
    #r.H <- matrix(0, ncol=p+1, nrow=H)
    weight_list <- seq(0, 1, length=H+2)[2:(H+1)]
    if(as.character(loss) == "ls"){
      A <- vector(mode = "list", length = length(qprob))
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
      # estimate bases
      r.H <- matrix(0, ncol=p+1, nrow=H)
      for (s in 1:length(qprob)){
        y.tilde <- rep(1, nrow(x))
        y.tilde[y < qy[s]] <- -1  #s
        A[[s]] <- M <- solve(n * cov.x.star/lambda+t(x.star) %*% x.star)
        C <- t(x.star) %*%  y.tilde
        r.H[s,] <- M %*% C
      }
    }
    if(as.character(loss)=="wls"){
      A <- vector(mode = "list", length = H)
      r.H <- matrix(0, ncol=p+1, nrow=H)
      if(sum(unique(y)) != 0)
        stop("response variable should be a binary type!")
      for (i in 1:H){
        W <- diag(c(ifelse(y==1, 1-weight_list[i], weight_list[i])))
        A[[i]] <- solve(n*cov.x.star/lambda + t(x.star) %*% W %*% x.star)
        C <- t(x.star) %*% W %*% y
        r.H[i,] <- A[[i]] %*% C
      }
    }
  }
  Working_mat <- t(r.H[,1:p])%*%r.H[,1:p]
  eigen.Mn <- eigen(Working_mat)

  newlist <- list("x" = x, "y" = y, "Mn"=Working_mat, "values" = eigen.Mn$values, "vectors" = eigen.Mn$vectors, "N"=n, "Xbar"=apply(x, 2, mean), "r"=r.H, "A"=A)
  class(newlist) <- "psdr"
  structure(class = "psdr", newlist)
  return(newlist)
}


#' @noRd
#' @export
print.psdr <- function(x, ...) {
  obj <- x
  #d <- list(x = obj$x, y = obj$y, Mn= obj$Mn, evalues = obj$values, evectors = obj$vectors)
  d <- list(Mn= obj$Mn, evalues = obj$values, evectors = obj$vectors)
  writeLines("psdr result:")
  print(d, ...)
  invisible(d)
}



#' @noRd
#' @export
plot.psdr <- function(x, dim=2, ...) {
  obj <- x
  if (!inherits(obj, "psdr"))
    stop("use only with \"psdr\" objects")
  temp <- obj$vectors
  obj_psdr <- obj$x %*% temp
  par(mfrow=c(ceiling(sqrt(dim)), ceiling(sqrt(dim))))
  par(mfrow=c(1,dim))
  for(d in 1:dim){
    plot(obj_psdr[,d], obj$y, type = "p", xlab = bquote(paste(hat(b)[.(d)]^T*X)), ylab  = expression(Y) , cex=.7,...)
    graphics::lines(lowess(obj_psdr[,d], obj$y), col="red", lwd=1)
  }
  par(mfrow=c(1,1))
}


