###################################
### Pre-embedded arbitrary loss ###
###################################

E <- function(object){UseMethod("E")}
fn_arbitrary <- function(object){UseMethod("fn_arbitrary")}
my.logistic <- function(object){UseMethod("my.logistic")}
my.hinge <- function(object){UseMethod("my.hinge")}
my.l2svm <- function(object){UseMethod("my.l2svm")}
my.lum <- function(object){UseMethod("my.lum")}
my.quantile <- function(object){UseMethod("my.quantile")}
my.asymLS <- function(object){UseMethod("my.asymLS")}
fn_arbitrary_loss <- function(object){UseMethod("fn_arbitrary_loss")}
fn_arbitrary_binary_loss <- function(object){UseMethod("fn_arbitrary_binary_loss")}
fn_arbitrary_nonlinear_loss <- function(object){UseMethod("fn_arbitrary_nonlinear_loss")}
fn_arbitrary_nonlinear_binary_loss <- function(object){UseMethod("fn_arbitrary_nonlinear_binary_loss")}
get.psi <- function(object){UseMethod("get.psi")}
psi.function <- function(object){UseMethod("psi.function")}
wvec <- function(object){UseMethod("wvec")}
kernel.function <- function(object){UseMethod("kernel.function")}
phix <- function(object){UseMethod("phix")}




E <- function (...) {eval(parse(text=paste(...,collapse=" ")))}


fn_arbitrary <- function(loss){
  E(loss)
}

my.logistic <- function(m,...){
  rslt <- log(1+exp(-m))
  return(rslt)
}

my.hinge <- function(m,...){
  rslt <- (1-m)*(as.numeric((1-m) > 0))
  return(rslt)
}

my.l2svm <- function(m,...){
  rslt <- ((1-m)*(as.numeric((1-m) > 0)))^2
  return(rslt)
}

my.lum <- function(m,a,c,...){
  rslt <- ((1-m)*(as.numeric(m < c/(1+c)))) + (1/(1+c))*({a/ ((1+c)*m-c+a) }^a)*(as.numeric(m >= c/(1+c)))
  return(rslt)
}

my.quantile <- function(u,prob,...){
  rslt <- abs(u)*prob*I(u>=0)+abs(u)*(1-prob)*I(u<=0)
  return(rslt)
}

my.asymLS <- function(u,prob,...){
  rslt <- (u)^2*prob*I(u>=0)+(u)^2*(1-prob)*I(u<=0)
  return(rslt)
}

############################
### Derivative functions ###
############################
##linear##

fn_arbitrary_loss <- function(x, y, theta, prob=0.5, lambda, loss,a=NULL,c=NULL){
  if(as.character(loss) == "my.quantile"){
    u <- y - x%*%theta
    ft <- fn_arbitrary(loss)
    losses <- ft(u,prob)
    loss.output <- lambda*mean(losses) + t(theta)%*%theta
  }else{
    m <- (x %*% theta)*y
    ft <- fn_arbitrary(loss)
    losses <- ft(m,a,c)
    loss.output <- lambda*mean(losses) + t(theta)%*%theta
  }
  return(loss.output)
}

fn_arbitrary_binary_loss <- function(x, y, prob=0.5, theta, lambda, loss,a=NULL,c=NULL){
  weight <- (1-prob)*(as.numeric(y==1)) + (prob)*(as.numeric(y==-1))
  if(as.character(loss)=="my.quantile"){
    u <- y - x%*%theta
    ft <- fn_arbitrary(loss)
    losses <- weight*ft(u,prob)
    loss.output <- lambda*mean(losses) + t(theta)%*%theta
  }else{
    m <- (x %*% theta)*y
    ft <- fn_arbitrary(loss)
    losses <- weight*ft(m,a,c)
    loss.output <- lambda*mean(losses) + t(theta)%*%theta
  }
  return(loss.output)
}

##nonlinear##
fn_arbitrary_nonlinear_loss <- function(x, y, theta, prob=NULL,lambda, loss,a=NULL,c=NULL){
  losses <- c()
  A <- t(x) %*% x
  if(as.character(loss)=="my.quantile"){
    u <- y - x%*%theta
    ft <- fn_arbitrary(loss)
    losses <- ft(u,prob)
    loss.output <- lambda*(sum(losses)/nrow(x)) + (t(theta) %*% A %*% theta)/nrow(x)
  }else{
    m <- (x %*% theta)*y
    ft <- fn_arbitrary(loss)
    losses <- ft(m,a,c)
    loss.output <- lambda*(sum(losses)/nrow(x)) + (t(theta) %*% A %*% theta)/nrow(x)
  }
  return(loss.output)
}


fn_arbitrary_nonlinear_binary_loss <- function(x, y, theta, prob,lambda, loss,a=NULL,c=NULL){
  losses <- c()
  A <- t(x) %*% x
  weight <- (1-prob)*(as.numeric(y==1)) + (prob)*(as.numeric(y==-1))
    m <- (x %*% theta)*y
    ft <- fn_arbitrary(loss)
    losses <- weight*ft(m,a,c)
    loss.output <- lambda*(sum(losses)/nrow(x)) + (t(theta) %*% A %*% theta)/nrow(x)
  return(loss.output)
}

###############################################
##auxilary functions for nonlinear prediction##
###############################################
get.psi <- function(x, y, k) {
    n <- nrow(x)
    x <- scale(x)
    bw <- 1/mean(as.numeric(dist(x)))^2 # bw parameter for kernel

    Kn <- svmpath::radial.kernel(x, x, bw)
    Qn <- diag(n) - matrix(1/n, n, n)

    eigen.psi <- eigen(Qn %*% Kn %*% Qn)
    Psi.new <- eigen.psi$vectors[,1:k, drop = F] # Psi
    l <- eigen.psi$values[1:k]             # eigen
    tmp.obj <- list("w"=Psi.new, "l"=l, "scaled.x"= x, "bw" = bw, "k" = k)
    class(tmp.obj) <- "npsdr"
    return(tmp.obj)
}


psi.function <- function(value, x, v, w, l, kernel.function, kernel.param){
  value <- matrix(value, 1, length(value))
  temp <- kernel.function(value, x, kernel.param)
  psi.value <- apply(w * c(temp - mean(temp)), 2, sum)/l
  psi.value %*% v
}

wvec <- function(x, y = y) {
    w <- NULL
    w[y == 1] <- 1-x
    w[y ==-1] <- x
    w
}

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


phix <- function(value, obj, order = 2) {
  psi.function <- psi.function
  x <- value
  v <- obj$evector
  w <- obj$obj.psi$w
  l <- obj$obj.psi$l
  kernel.function <- kernel.function
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

d2 <- function(Bhat, B) {
  # This function reports performance in terms of the Frobenius norm of the projection matrix difference

  return(norm(B%*%solve(t(B)%*%B)%*%t(B)-Bhat%*%solve(t(Bhat)%*%Bhat)%*%t(Bhat),"f"))
}
