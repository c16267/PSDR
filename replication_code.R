#
# This script contains code of all examples in the manuscript
# "upsdr: Unified Principal Sufficient Dimension Reduction R package"
# by Jungmin Shin, Artemiou Andreas, Seung Jun Shin
#


#=============================#
#### Load "psvmSDR" package ####
#=============================#
#install.packages("psvmSDR")
library("psvmSDR")

#=================================================#
#### Chapter 4 Basic Implementation           ####
#=================================================#

#=================================================#
#### 4.1 Functions for class 'psdr'            ####
#=================================================#

#===============================#
#### psdr() with loss='svm' #####
#===============================#
set.seed(100)
n <- 200;
p <- 5;
x <- matrix(rnorm(n*p, 0, 1), n, p)
y <-  x[,1] / (0.5 + (x[,2] + 1)^2) + 0.2*rnorm(n)
obj <- psdr(x, y)

#=====================#
#### print.psdr() #####
#=====================#
print(obj)

#=====================#
#### plot.psdr() #####
#=====================#
plot(obj)

#===============================#
#### psdr() with loss='wsvm' #####
#===============================#
y.binary <- sign(y)
obj_wsvm <- psdr(x, y.binary, loss="wsvm")
print(obj_wsvm)

#===============================================#
#### 3d plot for viewing the result of 'wsvm' #####
#==============================================#
plot(obj_wsvm)


#===============================#
#### psdr() with loss='mylogit' #####
#===============================#
mylogistic <- function(u, type="m") log(1+exp(-u))
obj_mylogistic <- psdr(x, y, loss="mylogistic")
print(obj_mylogistic)


#=====================#
#### crBIC()      #####
#=====================#
d.hat <- psdr_bic(obj)
print(d.hat)


#=================================================#
#### 4.3 Functions for class 'npsdr'            ####
#=================================================#
set.seed(1)
n <- 200;
p <- 5;
x <- matrix(rnorm(n*p, 0, 1), n, p)
y <-  log(x[,1]^2) * (x[,1]^2+x[,2]) + 0.2*rnorm(n)

#=====================#
#### linear case    #####
#=====================#
obj_lin <- psdr(x, y, loss="lssvm")
plot(obj_lin, main="Linear PLSSVM")

#=====================#
#### npsdr()      #####
#=====================#
obj_kernel <- npsdr(x, y, loss='lssvm', plot=FALSE)
print(obj_kernel)
plot(obj_kernel,  main="Nonlinear PLSSVM")


#=====================#
#### npsdrx()      #####
#=====================#
set.seed(1)
new.x <- matrix(rnorm(n*p, 0, 1), n, p)
new.y<-  log(new.x[,1]^2) * (new.x[,1]^2+new.x[,2]) + 0.2*rnorm(n)
reduced_data <- npsdr_x(object=obj_kernel, newdata=new.x, d = 2)

#==================================================================#
#### Linear Regression under the estimated central subspace    #####
#=================================================================#
#fit linear regression under the reduced dimension via npsdr
lm_npsdr <- lm(new.y ~ reduced_data[,c(1,2)])

#scatter plot y vs fitted value
par(mar=c(5,5,5,5), oma=c(1,1,1,1))
#par(mar=c(5.1, 4.1, 4.1, 2.1))
plot(scale(new.y), scale(lm_npsdr$fitted.values), pch=16, ylab=expression(paste("new", hat(Y))), xlab=expression(newY), cex.lab=1.5,
     ylim=c(-3,3), xlim=c(-3,3))
abline(a=0,b=1, col='red', lty=2, lwd=2)
legend('topright',lty=c(2),bg='white', col='red', lwd=2,legend=expression(y==x), cex=1.3)
grid(nx = NULL, ny = NULL, lty = 1, col = "gray", lwd = 1)


#=================================================#
#### 4.4 Functions for class 'rtpsdr'            ####
#=================================================#
n.sim <- 100
p <- 5
m <- 500 # batch size
N <- 10  # number of batches
norm_mat <- matrix(NA, ncol=4, nrow=n.sim)
colnames(norm_mat) <- c("rt_M", "psdr_M", "rt_r", "psdr_r")

for(j in 1:n.sim){
obj <- NULL
X <- matrix(NA, nrow=m*N, ncol=p)
Y <- c()

for (iter in 1:N){
  set.seed(iter+j)
  x <- matrix(rnorm(m*p, 0, 1), m, p)
  y <-  x[,1]/(0.5 + (x[,2] + 1)^2) + 0.2 * rnorm(m)
  #cumulative data
  X[(m*(iter-1)+1) : (iter*m), ] <- x
  Y <- c(Y,y)
  #SDR
  obj <- rtpsdr(x = x, y = y, obj=obj)
}
sdr_at_once <- psdr(x=X, y=Y, loss="lssvm")
print(paste("n.sim:", j))
}

#compare estimiation result
round(obj$evectors - sdr_at_once$evectors, 5)

#================================================#
#####real time weighted  least squares svm    #####
#================================================#
#load wisconsin breast cancer data from mclust package
#install.packages("mclust", force="TRUE")
library(mclust)
data("wdbc")
wisc <- wdbc
wisc <- wisc[1:550,]

n.sim <- 100
p <- ncol(wisc)-2
m <- 110# batch size
N <- 5  # number of batches
norm_mat <- matrix(NA, ncol=4, nrow=n.sim)
colnames(norm_mat) <- c("rt_M", "psdr_M", "rt_r", "psdr_r")


for(j in 1:n.sim){
  obj <- NULL
  X <- matrix(NA, nrow=m*N, ncol=p)
  Y <- c()

  for (iter in 1:N){
    set.seed(iter+j)
    #ind <- sample(1:550, m)
    ind <- (m*(iter-1)+1) : (iter*m)
    x.wisc <- matrix(unlist(wisc[ind, -c(1,2)]), ncol = 30)
    y.wisc <- (2*as.numeric(as.factor(unlist(wisc[ind, 2]))) - 3)

    #cumulative data
    X[(m*(iter-1)+1) : (iter*m), ] <- x.wisc
    Y <- c(Y,y.wisc)
    #SDR
    obj <- rtpsdr(x = x.wisc, y = y.wisc, obj=obj)
  }
  sdr_at_once <- psdr(x=X, y=Y, loss="wlssvm")

  norm_mat[j,1] <- norm(as.matrix(obj$M), "F")
  norm_mat[j,2] <- norm(as.matrix(sdr_at_once$M), "F")

  norm_mat[j,3] <- norm(obj$r, "F")
  norm_mat[j,4] <- norm(sdr_at_once$r, "F")

  print(paste("n.sim:", j))
}
round(obj$evectors[1:5, 1:5], 5)
round(sdr_at_once$evectors[1:5, 1:5], 5)

round(obj$r[1:5, 1:5], 5)
round(sdr_at_once$r[1:5, 1:5], 5)


round(apply(norm_mat, 2, "mean"), 5)
round(apply(norm_mat, 2, "sd")/sqrt(n.sim), 5)


#realtime
x.lsvm <- x.wisc %*% obj$evectors
par(mar=c(4,5,4,4), oma=c(1,2,1,1))
plot(x.lsvm[,1], x.lsvm[,2], type = "n", xlab = "Sufficient predictor 1", ylab  = "Sufficient predictor 2",
     main="Realtime PWLSSVM")
points(x.lsvm[y.wisc == 1,1], x.lsvm[y.wisc == 1,2], col = 2, pch = 16) #Malignant
points(x.lsvm[y.wisc != 1,1], x.lsvm[y.wisc != 1,2], col = 4, pch=16) #Benign
legend("top", legend=c("Benign", "Malignant"), col=c(4, 2), pch=c(16,16), cex=1.3)


#cumulative
x.lsvm <- X %*% sdr_at_once$evectors
par(mar=c(4,5,4,4), oma=c(1,2,1,1))
plot(x.lsvm[,1], x.lsvm[,2], type = "n", xlab = "Sufficient predictor 1", ylab  = "Sufficient predictor 2",
     main="Cumulative PWLSSVM")
points(x.lsvm[Y == 1,1], x.lsvm[Y == 1,2], col = 2, pch = 16) #Malignant
points(x.lsvm[Y != 1,1], x.lsvm[Y != 1,2], col = 4, pch=16) #Benign
legend("top", legend=c("Benign", "Malignant"), col=c(4, 2), pch=c(16,16), cex=1.3)



#============================================================#
##### Additional: Wisconsin Diagnostic Breast Cancer data #####
#=============================================================#

#=============================================#
##### linear weighted PLR        #####
#=============================================#
#load wisc dataset from data folder.

x.wisc <- matrix(unlist(wisc[,-c(1,2)]), ncol = 30)
y.wisc <- 2*as.numeric(as.factor(unlist(wisc[,2]))) - 3

obj_wisc <- psdr(x.wisc, y.wisc, loss='wlogit', h=20, lambda=0.1)
x.lsvm <- x.wisc %*% obj_wisc$evectors
par(mar=c(4,5,4,4), oma=c(1,2,1,1))
plot(x.lsvm[,1], x.lsvm[,2], type = "n", xlab = "Sufficient predictor 1", ylab  = "Sufficient predictor 2")
points(x.lsvm[y.wisc == 1,1], x.lsvm[y.wisc == 1,2], col = 2, pch = 16) #Malignant
points(x.lsvm[y.wisc != 1,1], x.lsvm[y.wisc != 1,2], col = 4, pch=16) #Benign
legend("top", legend=c("Benign", "Malignant"), col=c(4, 2), pch=c(16,16), cex=1.5)




#=============================================#
##### nonlinear weighted PLR        #####
#=============================================#
#=============================================#
# auxiliary functions 
#=============================================#
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


psi.function <- function(value, x, v, w, l, kernel.function, kernel.param){
  value <- matrix(value, 1, length(value))
  temp <- kernel.function(value, x, kernel.param)
  psi.value <- apply(w * c(temp - mean(temp)), 2, sum)/l
  rslt <- psi.value %*% v
  class(rslt) <- "npsdr"
  return(rslt)
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
nonlin_wisc_obj <- npsdr(x.wisc, y.wisc, loss="wsvm", h=10, lambda=1, max.iter = 30, plot=FALSE)
x.nsvm <- phix(x.wisc, nonlin_wisc_obj, d=2)  #recover to the original space
boxplot.default(x.nsvm[y.wisc == 1,1], x.nsvm[y.wisc != 1,1], xlab = "Y", axes = FALSE,
                ylab = expression(hat(phi)[1](x)))
axis(1, seq(0.5, 2.5, by = 0.5), c(NA, "+1", NA, "-1", NA)); axis(2, las = 1)



