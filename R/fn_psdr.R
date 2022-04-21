#'Gradient descent based sufficient dimension reduction machine
#'@description
#'A method for reducing data dimension using gradient descent optimization method
#'
#'@param x data matrix
#'@param y either continuous or (+1,-1) typed binary response vector
#'@param init.theta initial parameter theta
#'@param H the number of slicing
#'@param lambda hyperparameter for the loss function
#'@param h very small interval for calculating numerical derivatives for a given arbitrary loss function
#'@param delta learning rate
#'@param eps threshold for stopping iteration with respect to the magnitude of derivative
#'@param max.iter maximum iteration number for the optimization process
#'@param loss pre-specified loss functions are "logistic", svm","l2svm","lwpsvm", and user-defined loss function object also can be used formed by inside double quotation mark
#'@param a hyperparameter for LUM type loss function
#'@param c hyperparameter for LUM type loss function
#'@param stochastic user's can implement batch stochastic gradient descent method in case of large number of observations. default is set to FALSE.
#'@return An estimated matrix B for a sufficient dimension reduction will be returned
#'@author Jungmin Shin, \email{jungminshin@korea.ac.kr}, Seungjun Shin, \email{sjshin@korea.ac.kr}
#'@seealso \code{\link{plot.psdr}}
#'@examples
#'\donttest{
#'set.seed(1234)
#'n <- 300; p <- 5; H <- 10; lambda <- 0.1; eps <- 1.0e-6;
#'max.iter <- 10; h <- 1.0e-5; delta <- 3*1.0e-1;
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
#'my.logistic <- function(m,...){
#'rslt <- log(1+exp(-m))
#'return(rslt)
#'}
#'GD.psdr.solver(x, y, init.theta, H,lambda, h, delta, eps, max.iter, loss="svm")
#'GD.psdr.solver(x, y, init.theta, H,lambda, h, delta, eps, max.iter, loss="LUM", a=1, c=10^8)
#'GD.psdr.solver(x, y.binary, init.theta, H, lambda=0.5, h, delta, eps, max.iter, loss="my.logistic")
#'
#'#wisconsin breast cancer data example
#'wisc <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data", sep = ",")
#'x.wisc <- matrix(unlist(wisc[,-c(1,2)]), ncol = 30)
#'y.wisc <- y.wisc.binary <- 2*as.numeric(as.factor(unlist(wisc[,2]))) - 3
#'n <- length(y.wisc.binary) #n=569
#'p <- ncol(x.wisc)  #p=30
#'init.theta <- rnorm(sd=1,n=p)
#'wisc.obj <- GD.psdr.solver(x.wisc, y.wisc.binary, init.theta, H=20,lambda=0.1, h=1.0e-6, delta=0.5,
#'                           eps=10^-4, max.iter=30, loss="wlogistic")
#'value.lsvm <- wisc.obj$values
#'lsvm <- round(wisc.obj$vectors,3)
#'x.lsvm <- x.wisc %*% lsvm
#'plot(x.lsvm[,1], x.lsvm[,2], type = "n", xlab = "1st predictor", ylab  = "2nd predictor")
#'points(x.lsvm[y.wisc.binary == 1,1], x.lsvm[y.wisc.binary == 1,2], col = 4, pch = "+")
#'points(x.lsvm[y.wisc.binary != 1,1], x.lsvm[y.wisc.binary != 1,2], col = 2)
#'}
#'@import stats MASS graphics
#'@export


##################################################
### main function ################################
##################################################

GD.psdr.solver <- function(x, y, init.theta, H, lambda, h, delta, eps, max.iter=NULL,
                           loss, a=NULL, c=NULL, stochastic=FALSE) {

  if(sum(as.character(loss) == c("ls", "wls")) == 0){
    if(!is.matrix(x) & !is.data.frame(x))
      stop("x must be a matrix or dataframe.")
    if(ncol(as.matrix(y)) != 1)
      stop("y must be a univariate.")
    if(length(init.theta) != ncol(x))
      stop("a dimension of the initial theta must match that of input x.")
    if(as.numeric(is.null(delta)) == 1 | delta < 0)
      stop("Delta is a learning rate which should be specified as a positve value.")
    if(as.character(loss) == "LUM" & sum(is.null(c(a,c))) >= 1)
      stop("A LUM loss fucntion needs both hyperparameter a and c.")
  }
  if(H <= 1)
    stop("H should be greater than 1.")
  if(lambda <= 0)
    stop("lambda should be a positive numeric.")

  if (length(y) != nrow(x))     #check lengths
    stop("The response and predictors have different number of observations")

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
  type.list <- c("svm","logistic","l2svm", "lwpsvm", "LUM", "quantile","wlogistic","l2wsvm","wLUM")
  type.list2 <- c("ls","wls")

  if(sum(as.character(loss) == type.list) != 0){
    w.init <- matrix(init.theta, nrow=p, ncol=length(qprob))
    w.final <- matrix(0, nrow=p, ncol=length(qprob))
    if(as.character(loss) == "svm"){
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

    if(as.character(loss) == "lwpsvm"){
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
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
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

    psi <- t(inv.sd.x) %*% w.final
    Mn <- matrix(0, p, p)
    for (h in 1:length(qprob)) Mn <- Mn + psi[,h, drop = F] %*% t(psi[,h, drop = F])
    eigen.Mn <- eigen(Mn)
    newlist <- list("Mn"=Mn, "values" = eigen.Mn$values, "vectors" = eigen.Mn$vectors)
    return(newlist)
    class(newlist) <- "psdr"
  }

  else if(sum(as.character(loss) == type.list) == 0 & sum(as.character(loss)==type.list2)==0){
    if(h > 1.0e-3)
      warning("An infinitesimal interval h should be smaller than 1.0e-3.")
    ft <- E(loss)
    grid.m <- seq(-2,2,length=100)
    plot(grid.m, ft(grid.m,a,c,prob=0.5), type="l", xlab="margin", ylab="loss")
    message("loss function must be a convex function")
    w.init <- matrix(init.theta, nrow=p, ncol=length(qprob))
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
    newlist <- list("Mn"=Mn, "values" = eigen.Mn$values, "vectors" = eigen.Mn$vectors)
    return(newlist)
    class(newlist) <- "psdr"
  }

  else if(sum(as.character(loss) == type.list2) != 0){
    r.H <- matrix(0, ncol=p+1, nrow=H)
    weight_list <- seq(0, 1, length=H+2)[2:(H+1)]
    if(as.character(loss) == "ls"){
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
      # estimate bases
      r.H <- matrix(0, ncol=p+1, nrow=H)
      for (s in 1:length(qprob)){
        y.tilde <- rep(1, nrow(x))
        y.tilde[y < qy[s]] <- -1  #s
        M <- solve(n * cov.x.star/lambda+t(x.star) %*% x.star)
        C <- t(x.star) %*%  y.tilde
        r.H[s,] <- M %*% C
      }
    }
    if(as.character(loss)=="wls"){
      if(sum(unique(y)) != 0)
        stop("response variable should be a binary type!")
      for (pi in weight_list){
        W <- diag(c(ifelse(y==1, 1-pi, pi)))
        M <- solve(n*cov.x.star/lambda + t(x.star) %*% W %*% x.star)
        C <- t(x.star) %*% W %*% y
        r.H[which(weight_list==pi),] <- M %*% C
      }
    }
  }
  eigen.Mn <- eigen(t(r.H[,1:p])%*%r.H[,1:p])
  newlist <- list("Mn"=r.H, "values" = eigen.Mn$values, "vectors" = eigen.Mn$vectors)
  return(newlist)
  class(newlist) <- "psdr"
}


