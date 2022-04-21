#'Gradient descent based nonlinear sufficient dimension reduction machine via kernel trick
#'@description
#'A method for reducing a dimension of the data sufficiently using gradient descent optimization method
#'
#'@param x data matrix
#'@param y either continuous or (+1,-1) typed binary response vector
#'@param H the number of slicing
#'@param h very small interval for calculating numerical derivatives for a given arbitrary loss function
#'@param lambda hyperparameter for the loss function
#'@param delta learning rate
#'@param k number of basis functions for a kernel trick
#'@param eps threshold for stopping iteration with respect to the magnitude of derivative
#'@param max.iter maximum iteration number for the optimization process
#'@param loss pre-specified loss functions are "logistic", svm","l2svm","lwpsvm", and user-defined loss function object also can be used formed by inside double quotation mark
#'@param a a hyperparameter for the LUM loss function
#'@param c a hyperparameter for the LUM loss function
#'@return An estimated matrix B for a sufficient dimension reduction will be returned
#'@author Jungmin Shin, \email{jungminshin@korea.ac.kr}, Seungjun Shin, \email{sjshin@korea.ac.kr}
#'@seealso \code{\link{pred.npsdr}}
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
#'my.hinge <- function(m,...){
#'  rslt <- (1-m)*(as.numeric((1-m) > 0))
#'  return(rslt)
#'}
#'GD.npsdr.solver(x, y, H, h, lambda, delta, k=floor(length(y)/3), eps,
#'                max.iter, loss="my.logistic")
#'GD.npsdr.solver(x, y.binary, H=20, h, lambda, delta, k=floor(length(y)/3), eps,
#'                max.iter, loss="svm") }
#'@import stats kernlab svmpath
#'@export


GD.npsdr.solver <- function(x, y, H, h, lambda, delta, k = floor(length(y)/3), eps, max.iter, loss, a=NULL, c=NULL)
{
  if(sum(as.character(loss) == c("ls", "wls")) == 0){
    if(!is.matrix(x) & !is.data.frame(x))
      stop("x must be a matrix or dataframe.")
    if(ncol(as.matrix(y)) != 1)
      stop("y must be a univariate.")
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


  psi.gen <- get.psi(x,y,k)
  Psi.new <- psi.gen$w   #n*k
  n <- nrow(Psi.new)
  p <- ncol(Psi.new)
  bar.x <- apply(Psi.new, 2, mean)
  x.star <- cbind(t(t(Psi.new)-bar.x), -1)
  cov.x.star <- cov(x.star)
  init.theta <- rnorm(sd=1,n=p)
  step <- 1/H
  pi.grid <- seq(step, 1-step, by = step)

  # generate y.tilde
  qprob <- (1:(H-1))/H
  qy <- stats::quantile(y, qprob)

  theta.new <- rep(0,p)
  w.init <- matrix(init.theta, nrow=p, ncol=length(qprob))
  w.final <- matrix(0, nrow=p, ncol=length(qprob))
  eigen.mat <- diag(1,p,p)

  type.list <- c("svm","logistic","l2svm", "nwpsvm", "LUM", "wlogistic","wl2svm","wLUM","quantile")
  type.list2 <- c("ls","wls")

  if(sum(as.character(loss) == type.list) != 0){
    if(as.character(loss) == "svm"){
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
      for (s in 1:length(qprob)) {
        y.tilde.new <- rep(1, nrow(Psi.new))
        y.tilde.new[y < qy[s]] <- -1  #s
      #  pos.rate <- sum(y.tilde.new==1)/nrow(Psi.new)
      #  neg.rate <- sum(y.tilde.new==-1)/nrow(Psi.new)
        for(iter in 1:max.iter){
          #if(stochastic == TRUE){  ###stratified sampling
          #  set.seed(iter)
          #  pos.ind <- sample(which(y.tilde.new==1), ceiling((nrow(Psi.new)/(log(nrow(Psi.new))))*pos.rate))
          #  neg.ind <- sample(which(y.tilde.new==-1), floor((nrow(Psi.new)/(log(nrow(Psi.new))))*neg.rate))
          #  ind <- c(pos.ind,neg.ind)
          #  Psi <- Psi.new[ind,]
          #  y.tilde <- y.tilde.new[ind]
          #}else{
            Psi <- Psi.new
            y.tilde <- y.tilde.new
          #}
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi
          for (k in 1:p){
            margin.v <- Psi %*% w[,s] * y.tilde #s
            deriv <- -Psi[,k]*y.tilde*as.numeric(I((1-margin.v)>0)) #k
            derivative.j <- lambda*mean(deriv) + 2*(1/nrow(Psi))*(A[k,]%*%w[,s])  ##k,s
            theta.new[k] <- w[k,s] -  delta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new
          w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
          if(max(abs(deriv)) < eps)
            break
        }
        w.final[,s] <- w[,s] #s,s
      }
    }
    if(as.character(loss) == "l2svm"){
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
      for (s in 1:length(qprob)) {
        y.tilde.new <- rep(1, nrow(Psi.new))
        y.tilde.new[y < qy[s]] <- -1  #s
      #  pos.rate <- sum(y.tilde.new==1)/nrow(Psi.new)
      #  neg.rate <- sum(y.tilde.new==-1)/nrow(Psi.new)
        for(iter in 1:max.iter){
          #if(stochastic == TRUE){  ###stratified sampling
          #  set.seed(iter)
          #  pos.ind <- sample(which(y.tilde.new==1), ceiling((nrow(Psi.new)/(log(nrow(Psi.new))))*pos.rate))
          #  neg.ind <- sample(which(y.tilde.new==-1), floor((nrow(Psi.new)/(log(nrow(Psi.new))))*neg.rate))
          #  ind <- c(pos.ind,neg.ind)
          #  Psi <- Psi.new[ind,]
          #  y.tilde <- y.tilde.new[ind]
          #}else{
            Psi <- Psi.new
            y.tilde <- y.tilde.new
          #}
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi
          for (k in 1:p){
            margin.v <- Psi %*% w[,s] * y.tilde #s
            deriv <- 2*(1-margin.v)*(-Psi[,k]*y.tilde)*as.numeric(I((1-margin.v)>0)) #k
            derivative.j <- lambda*mean(deriv) + 2*(1/nrow(Psi))*A[k,]%*%w[,s]  ##k,s
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
    if(as.character(loss) == "wl2svm"){
      if(sum(unique(y)) != 0)
        stop("response variable should be a binary type!")
      y.new <- y
      for (s in 1:length(qprob)) {
      #  pos.rate <- sum(y.new==1)/nrow(Psi.new)
      #  neg.rate <- sum(y.new==-1)/nrow(Psi.new)
        for(iter in 1:max.iter){
          #if(stochastic == TRUE){  ###stratified sampling
          #  set.seed(iter)
          #  pos.ind <- sample(which(y.new==1), ceiling((nrow(Psi.new)/(log(nrow(Psi.new))))*pos.rate))
          #  neg.ind <- sample(which(y.new==-1), floor((nrow(Psi.new)/(log(nrow(Psi.new))))*neg.rate))
          #  ind <- c(pos.ind,neg.ind)
          #  Psi <- Psi.new[ind,]
          #  y.bi <- y.new[ind]
          #}else{
            Psi <- Psi.new
            y.bi <- y.new
          #}
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi

          for (k in 1:p){
            margin.v <- Psi %*% w[,s] * y.bi #s
            weight <- (1-pi.grid[s])*as.numeric(y.bi==1)+(pi.grid[s])*(as.numeric(y.bi==-1)) #s,s
            deriv <- 2*(1-margin.v)*(-Psi[,k]*y.bi)*as.numeric(I((1-margin.v)>0)) #k
            derivative.j <- lambda*mean(deriv) + 2*(1/nrow(Psi))*A[k,]%*%w[,s]  ##k,s
            theta.new[k] <- w[k,s] -  delta*derivative.j  ##k,k,s
          }

          w[,s] <- theta.new  #c
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
        y.tilde.new <- rep(1, nrow(Psi.new))
        y.tilde.new[y < qy[s]] <- -1  #s
        #pos.rate <- sum(y.tilde.new==1)/nrow(Psi.new)
        #neg.rate <- sum(y.tilde.new==-1)/nrow(Psi.new)
        for(iter in 1:max.iter){
        #  if(stochastic == TRUE){  ###stratified sampling
        #    set.seed(iter)
        #    pos.ind <- sample(which(y.tilde.new==1), ceiling((nrow(Psi.new)/(log(nrow(Psi.new))))*pos.rate))
        #    neg.ind <- sample(which(y.tilde.new==-1), floor((nrow(Psi.new)/(log(nrow(Psi.new))))*neg.rate))
        #    ind <- c(pos.ind,neg.ind)
        #    Psi <- Psi.new[ind,]
        #    y.tilde <- y.tilde.new[ind]
        #  }else{
            Psi <- Psi.new
            y.tilde <- y.tilde.new
         # }
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi

          for (k in 1:p){
            margin.v <- Psi %*% w[,s] * y.tilde #s
            deriv <- -Psi[,k]*y.tilde*(1/(1+exp(margin.v)))#k
            derivative.j <- lambda*mean(deriv) + 2*(1/nrow(Psi))*A[k,]%*%w[,s]  ##k,s
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
      if(sum(unique(y)) != 0)
        stop("response variable should be a binary type!")
      y.new <- y
      for (s in 1:length(qprob)){
        #pos.rate <- sum(y.new==1)/nrow(Psi.new)
        #neg.rate <- sum(y.new==-1)/nrow(Psi.new)
        for(iter in 1:max.iter){
        #  if(stochastic == TRUE){  ###stratified sampling
        #    set.seed(iter)
        #    pos.ind <- sample(which(y.new==1), ceiling((nrow(Psi.new)/(log(nrow(Psi.new))))*pos.rate))
        #    neg.ind <- sample(which(y.new==-1), floor((nrow(Psi.new)/(log(nrow(Psi.new))))*neg.rate))
        #    ind <- c(pos.ind,neg.ind)
        #    Psi <- Psi.new[ind,]
        #    y.bi <- y.new[ind]
        #  }else{
            Psi <- Psi.new
            y.bi <- y.new
        # }
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi
          for (k in 1:p){
            margin.v <- Psi %*% w[,s] * y.bi #s
            weight <- (1-pi.grid[s])*as.numeric(y.bi==1)+(pi.grid[s])*(as.numeric(y.bi==-1)) #s,s
            deriv <- weight*(-Psi[,k])*y.bi*(1/(1+exp(margin.v)))   #k
            derivative.j <- lambda*mean(deriv) + 2*(1/nrow(Psi))*A[k,]%*%w[,s]  ##k,s
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

    if(as.character(loss) == "nwpsvm"){
      if(sum(unique(y)) != 0)
        stop("response variable should be a binary type!")
      y.new <- y
      for (s in 1:length(qprob)) {
        #pos.rate <- sum(y.new==1)/nrow(Psi.new)
        #neg.rate <- sum(y.new==-1)/nrow(Psi.new)
        for(iter in 1:max.iter){
        #  if(stochastic == TRUE){  ###stratified sampling
        #    set.seed(iter)
        #    pos.ind <- sample(which(y.new==1), ceiling((nrow(Psi.new)/(log(nrow(Psi.new))))*pos.rate))
        #    neg.ind <- sample(which(y.new==-1), floor((nrow(Psi.new)/(log(nrow(Psi.new))))*neg.rate))
        #    ind <- c(pos.ind,neg.ind)
        #    Psi <- Psi.new[ind,]
        #    y.bi <- y.new[ind]
        #  }else{
            Psi <- Psi.new
            y.bi <- y.new
        #  }
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi

          for (k in 1:p){
            margin.v <- Psi %*% w[,s] * y.bi #s
            weight <- (1-pi.grid[s])*as.numeric(y.bi==1)+(pi.grid[s])*(as.numeric(y.bi==-1)) #s,s
            deriv <- -weight*Psi[,k]*y.bi*as.numeric(I((1-margin.v)>0)) #k
            derivative.j <- lambda*mean(deriv) + 2*(1/nrow(Psi))*A[k,]%*%w[,s]  ##k,s
            theta.new[k] <- w[k,s] -  delta*derivative.j  ##k,k,s
          }

          w[,s] <- theta.new  #c
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
        y.tilde.new <- rep(1, nrow(Psi.new))
        y.tilde.new[y < qy[s]] <- -1  #s
        #pos.rate <- sum(y.tilde.new==1)/nrow(Psi.new)
        #neg.rate <- sum(y.tilde.new==-1)/nrow(Psi.new)
        for(iter in 1:max.iter){
        #  if(stochastic == TRUE){  ###stratified sampling
        #    set.seed(iter)
        #    pos.ind <- sample(which(y.tilde.new==1), ceiling((nrow(Psi.new)/(log(nrow(Psi.new))))*pos.rate))
        #    neg.ind <- sample(which(y.tilde.new==-1), floor((nrow(Psi.new)/(log(nrow(Psi.new))))*neg.rate))
        #    ind <- c(pos.ind,neg.ind)
        #    Psi <- Psi.new[ind,]
        #    y.tilde <- y.tilde.new[ind]
        #  }else{
            Psi <- Psi.new
            y.tilde <- y.tilde.new
        #  }
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi

          for (k in 1:p){
            margin.v <- Psi %*% w[,s] * y.tilde #s
            deriv <- -Psi[,k]*y.tilde*as.numeric(margin.v < c/(1+c)) +
              (a/(1+c))*{(a/((1+c)*margin.v-c+a))^(a-1)}*( (1+c)*a*(-Psi[,k]*y.tilde)) / {(1+c)*margin.v-c+a}^2*
              as.numeric(margin.v >= c/(1+c))
            derivative.j <- lambda*mean(deriv) + 2*(1/nrow(Psi))*A[k,]%*%w[,s]  ##k,s
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
      for (s in 1:length(qprob)){
        #pos.rate <- sum(y.new==1)/nrow(Psi.new)
        #neg.rate <- sum(y.new==-1)/nrow(Psi.new)
        for(iter in 1:max.iter){
        # if(stochastic == TRUE){  ###stratified sampling
        #    set.seed(iter)
        #    pos.ind <- sample(which(y.new==1), ceiling((nrow(Psi.new)/(log(nrow(Psi.new))))*pos.rate))
        #    neg.ind <- sample(which(y.new==-1), floor((nrow(Psi.new)/(log(nrow(Psi.new))))*neg.rate))
        #    ind <- c(pos.ind,neg.ind)
        #    Psi <- Psi.new[ind,]
        #    y.bi <- y.new[ind]
        #  }else{
            Psi <- Psi.new
            y.bi <- y.new
        #  }
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi
          for (k in 1:p){
            margin.v <- Psi %*% w[,s] * y.bi #s
            weight <- (1-pi.grid[s])*as.numeric(y.bi==1)+(pi.grid[s])*(as.numeric(y.bi==-1)) #s,s
            deriv <- {-Psi[,k]*y.bi*as.numeric(margin.v < c/(1+c)) +
              (a/(1+c))*{(a/((1+c)*margin.v-c+a))^(a-1)}*( (1+c)*a*(-Psi[,k]*y.bi)) / {(1+c)*margin.v-c+a}^2*
              as.numeric(margin.v >= c/(1+c))}*weight
            derivative.j <- lambda*mean(deriv) + 2*(1/nrow(Psi))*A[k,]%*%w[,s]  ##k,s
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

    if(as.character(loss) == "quantile"){
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
      for (s in 1:length(pi.grid)) {
        for(iter in 1:max.iter){
        #  if(stochastic == TRUE){  ###stratified sampling
        #    set.seed(iter)
        #    ind <- sample(length(y), sqrt(p*length(y))) #sqrt(np)
        #    Psi <- Psi.new[ind,]
        #    y.new <- y[ind]
        #  }else{
            Psi <- Psi.new
            y.new <- y
        #  }
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi
          for (k in 1:p){
            u <- y.new - Psi%*%w[,s]   #s
            #derivative.j <- 2*w[k,s]+lambda*sum((-z[,k]*pi.grid[s])*as.numeric(I(u>0))+(z[,k]*(1-pi.grid[s]))*as.numeric(I(u<0)))*(1/length(y.new)) #k,s,k,s,k,s
            derivative.j <- 2*(1/nrow(Psi))*w[k,s]+lambda*(1/length(y.new))*sum(-Psi[,k]*{pi.grid[s]*as.numeric(I(u>0))+(1-pi.grid[s])*as.numeric(I(u<=0))})
            theta.new[k] <- w[k,s] -  delta*derivative.j  ##k,k,s
          }
          #print(derivative.j)
          w[,s] <- theta.new  #s
          w.init <- matrix(theta.new, nrow=p, ncol = length(pi.grid))
        }
        w.final[,s] <- w[,s] #s,s
      }
    }

    Mn <- matrix(0, p, p)
    for (h in 1:length(qprob)) Mn <- Mn + w.final[,h, drop = F] %*% t(w.final[,h, drop = F])
    result <- eigen(Mn)
    v <- result$vectors
    u <- result$values
    obj <- list(evector = v, evalue = u, obj.psi = psi.gen)
    class(obj) <- "npsdr"
    return(obj)
  }

  else if(sum(as.character(loss) == type.list) == 0 & sum(as.character(loss)==type.list2)==0){
    if(h > 1.0e-3)
      warning("An infinitesimal interval h should be smaller than 1.0e-3.")
    ft <- E(loss)
    grid.m <- seq(-2,2,length=100)
    plot(grid.m, ft(grid.m,a,c), type="l", xlab="margin", ylab="loss")
    message("loss function must be a convex function")

    ## arbitrary loss continuous##

    if(sum(unique(y)) != 0){
      ## main for loop ##
      for(s in 1:length(qprob)){
        y.tilde.new <- rep(1, nrow(Psi.new))
        y.tilde.new[y < qy[s]] <- -1  #s
        #pos.rate <- sum(y.tilde.new==1)/nrow(Psi.new)
        #neg.rate <- sum(y.tilde.new==-1)/nrow(Psi.new)
        w <- w.init
        for(iter in 1:max.iter){
         # if(stochastic == TRUE){  ###stratified sampling
         #    set.seed(iter)
         #    pos.ind <- sample(which(y.tilde.new==1), ceiling((nrow(Psi.new)/(log(nrow(Psi.new))))*pos.rate))
          #  neg.ind <- sample(which(y.tilde.new==-1), floor((nrow(Psi.new)/(log(nrow(Psi.new))))*neg.rate))
          #  ind <- c(pos.ind,neg.ind)
          #  Psi <- Psi.new[ind,]
          #  y.tilde <- y.tilde.new[ind]
          #}else{
            Psi <- Psi.new
            y.tilde <- y.tilde.new
          #}
          n <- nrow(Psi)
          derivative.vec <- rep(0,p)
          for (k in 1:p){
            star <- (fn_arbitrary_nonlinear_loss(Psi, y.tilde,theta=w[,s]+h*eigen.mat[,k],lambda,loss,a,c,prob=pi.grid[s])-
                       fn_arbitrary_nonlinear_loss(Psi, y.tilde,theta=w[,s]-h*eigen.mat[,k],lambda,loss,a,c,prob=pi.grid[s]))
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
      for(s in 1:length(qprob)){
        #pos.rate <- sum(y.new==1)/nrow(Psi.new)
        #neg.rate <- sum(y.new==-1)/nrow(Psi.new)
        w <- w.init
        for(iter in 1:max.iter){
        #  if(stochastic == TRUE){  ###stratified sampling
        #    set.seed(iter)
        #    pos.ind <- sample(which(y.new==1), ceiling((nrow(Psi.new)/(log(nrow(Psi.new))))*pos.rate))
        #    neg.ind <- sample(which(y.new==-1), floor((nrow(Psi.new)/(log(nrow(Psi.new))))*neg.rate))
        #    ind <- c(pos.ind,neg.ind)
        #    Psi <- Psi.new[ind,]
        #    y.bi <- y.new[ind]
        #  }else{
            Psi <- Psi.new
            y.bi <- y.new
        #  }
          n <- nrow(Psi)
          derivative.vec <- rep(0,p)
          for (k in 1:p){
            derivative.vec[k] <- (fn_arbitrary_nonlinear_binary_loss(Psi, y.bi, prob=pi.grid[s],theta=w[,s]+h*eigen.mat[,k],lambda,loss,a,c)-
                                    fn_arbitrary_nonlinear_binary_loss(Psi, y.bi ,prob=pi.grid[s],theta=w[,s]-h*eigen.mat[,k],lambda,loss,a,c)) / (2*h)
            theta.new[k] <- w[k,s] - delta * derivative.vec[k]
          }
          w[,s] <- theta.new
          if(mean(abs(derivative.vec)) < eps)
          break
        }
        w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
        w.final[,s] <- w[,s]
      }
    }

    Mn <- matrix(0, p, p)
    for (h in 1:length(qprob)) Mn <- Mn + w.final[,h, drop = F] %*% t(w.final[,h, drop = F])
    result <- eigen(Mn)
    v <- result$vectors
    u <- result$values
    obj <- list(evector = v, evalue = u, obj.psi = psi.gen)
    class(obj) <- "npsdr"
    return(obj)
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
        y.tilde <- rep(1, nrow(Psi.new))
        y.tilde[y < qy[s]] <- -1  #s
        M <- solve(n * cov.x.star/lambda+t(x.star) %*% x.star)
        C <- t(x.star) %*%  y.tilde
        r.H[s,] <- M %*% C
      }
    }
    if(as.character(loss)=="wls"){
      if(sum(unique(y)) != 0)
        stop("response variable should be a binary type!")
      #y.binary <- y
      for (pi in weight_list){
        W <- diag(c(ifelse(y==1, 1-pi, pi)))
        M <- solve(n*cov.x.star/lambda + t(x.star) %*% W %*% x.star)
        C <- t(x.star) %*% W %*% y
        r.H[which(weight_list==pi),] <- M %*% C
      }
    }
    result <- eigen(t(r.H[,1:p])%*%r.H[,1:p])
    v <- result$vectors
    u <- result$values
    newlist <- list(evector=v, evalue=u, obj.psi = psi.gen)
    #newlist <- list("Mn"=r.H, "values" = eigen.Mn$values, "vectors" = eigen.Mn$vectors)
    return(newlist)
    class(newlist) <- "npsdr"
  }
}




