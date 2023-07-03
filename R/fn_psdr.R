psdr <- function(x, y, init=NULL, H=NULL, lambda=NULL, delta=NULL, h=1.0e-5, eps=1.0e-5, max.iter=NULL,
                 loss=NULL, a=NULL, c=NULL, stochastic=FALSE) {
  if(sum(as.character(loss) == c("ls", "wls")) == 0){
    if(!is.matrix(x) & !is.data.frame(x))
      stop("x must be a matrix or dataframe.")
    if(ncol(as.matrix(y)) != 1)
      stop("y must be a univariate.")
    if(is.null(delta) == T ){
      warning("Delta is a learning rate which should be specified as a positve value. Delta is set to 0.1 as a default")
      delta <- 0.1
    }

    if(is.null(loss)==T){
      warning("Loss function is set to hinge loss as a default.")
      loss <- 'svm'
    }

    if(as.character(loss) == "LUM" ){
      if(is.null(a) == T)
        warning("A LUM loss fucntion needs both hyperparameter a and c. Default values are applied; a=1, c=10^8")
      a <- 1; c<- 10^8
    }

    if(is.null(H) == T){
      warning("H is set to 10 as a default.")
      H <- 10
    }
    if(is.null(h) == T){
      warning("h is set to 1.0e-5 as a default.")
      h <- 1.0e-5
    }
    if(is.null(max.iter) == T){
      warning("max.iter is set to 30 as a default.")
      max.iter <- 30
    }
    if(is.null(init) == T){
      init <- rnorm(sd=1,n=p)
      warning("initial parameter is generated from N(0,1)")
    }
    if(length(init) != ncol(x))
      stop("a dimension of the initial theta must match that of input x.")
  }


  if(is.null(lambda) == T){
    warning("lambda is set to 0.1 as a default.")
    lambda <- 0.1
  }
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
  type.list <- c("svm","logistic","l2svm", "lwpsvm", "LUM", "quantile","asymls", "wlogistic","l2wsvm","wLUM")
  type.list2 <- c("ls","wls")

  if(sum(as.character(loss) == type.list) != 0){
    w.init <- matrix(init, nrow=p, ncol=length(qprob))
    w.final <- matrix(0, nrow=p, ncol=length(qprob))

    if(as.character(loss) == "svm"){
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")

      for (s in 1:length(qprob)) {
        y.tilde.new <- rep(1, nrow(x))
        y.tilde.new[y < qy[s]] <- -1  #s
        pos.rate <- sum(y.tilde.new==1)/nrow(z.new)
        neg.rate <- sum(y.tilde.new==-1)/nrow(z.new)

        if(stochastic == TRUE ){
          for(iter in 1:max.iter){
            set.seed(iter)
            m <- nrow(x)
            mini_batches <- list()
            rand_sample <- c(sample(m))
            shuffled_X <- matrix(z.new[rand_sample,], ncol=p)
            shuffled_Y <- matrix(y.tilde.new[rand_sample], ncol=1)

            pos.ind <- sample(which(y.tilde.new==1), ceiling((nrow(z.new)/(log(nrow(z.new))))*pos.rate))
            neg.ind <- sample(which(y.tilde.new==-1), floor((nrow(z.new)/(log(nrow(z.new))))*neg.rate))
            ind <- c(pos.ind,neg.ind)
            z <- shuffled_X[ind,]
            y.tilde <- shuffled_Y[ind]
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
        }else{
          for(iter in 1:max.iter){
            z <- z.new
            y.tilde <- y.tilde.new
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
      warning("An infinitesimal interval h should be smaller than 1.0e-3.")
    ft <- E(loss)
    grid.m <- seq(-2,2,length=100)
    plot(grid.m, ft(grid.m,a,c,prob=0.5), type="l", xlab="margin", ylab="loss")
    message("loss function must be a convex function")
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




