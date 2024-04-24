#'Pre-embedded arbitrary loss
#' @noRd
#'
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
wvec <- function(object){UseMethod("wvec")}
d2 <- function(object){UseMethod("d2")}
E <- function(object){UseMethod("E")}


# get.psi <- function(object){UseMethod("get.psi")}
# psi.function <- function(object){UseMethod("psi.function")}
# kernel.function <- function(object){UseMethod("kernel.function")}
# phix <- function(object){UseMethod("phix")}


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


wvec <- function(x, y = y) {
  w <- NULL
  w[y == 1] <- 1-x
  w[y ==-1] <- x
  w
}

d2 <- function(Bhat, B) {
  # This function reports performance in terms of the Frobenius norm of the projection matrix difference

  return(norm(B%*%solve(t(B)%*%B)%*%t(B)-Bhat%*%solve(t(Bhat)%*%Bhat)%*%t(Bhat),"f"))
}






