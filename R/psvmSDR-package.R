#' A unified package for a principal sufficient dimension reduction method
#'
#' Principal Sufficient Dimension Reduction methods
#'
#' This package estimates the basis of the central subspace for both regression and classification problems
#' with underlying assumption of linear or nonlinear relationship between predictors and response variable.
#' This package implements a various principal sufficient dimension reduction methods (often called 'Principal Machines (PM)') in a unified fashion.
#' It ranges from linear, nonlinear dimension reduction to real-time cases where the data is collected sequentially.
#'
#' \tabular{ll}{ Package: \tab psvmSDR\cr Type: \tab Package\cr Version: \tab
#' 1.0.0\cr Date: \tab 2024-05-01\cr License: \tab GPL-2 \cr }
#' Very simple to use. Accepts \code{x,y} data for regression (or classification) models, and
#' produces the basis of the central subspace, which has a lower rank to the original data matrix.
#' The main 3 functions are: \code{psdr} for a linear principal machines (PM), \cr
#' \code{npsdr} is for a nonlinear PM \cr \code{rtpsdr} is for a real-timeprincipal least square SVM;
#'
#' @name psvmSDR-package
#' @docType package
#' @author Jungmin Shin, Seung Jun Shin, Andreas Artemiou \cr Maintainer:
#' Jungmin Shin \email{jungminshin@korea.ac.kr}
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
#'
#' @import graphics svmpath
#' @aliases psvmSDR-package
#' @srrVerbose FALSE
